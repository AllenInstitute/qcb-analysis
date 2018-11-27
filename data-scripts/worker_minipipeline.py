import os
import json
import errno
import pickle
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio

if __name__ == "__main__":

    #
    # External parameters
    #

    parser = argparse.ArgumentParser(description="Interface for downloading/exploring data from database and extracting features for QCB")
    parser.add_argument("-m", "--mode", help="feature (feature extraction) \
                                            contrast (adjust the contrast of images of a particular dataset)\
                                            config (generate the config json file based on a csv table)", required=True)
    parser.add_argument("-c", "--config", nargs="?", type=str, default="", const="", required=False)
    parser.add_argument("-d", "--dataset", nargs="?", type=str, default="", const="", required=False)
    parser.add_argument("-s", "--start", nargs="?", type=int, default=0, const=0, required=False)
    args = vars(parser.parse_args())

    #
    # Auxiliar functions
    #

    def get_orthogonal_projs(img):

        nz = img.shape[0]
        img_proj_xy = img.max(axis=0)
        img_proj_xz = img.max(axis=1)
        img_proj_yz = img.max(axis=2)
        img_offset  = np.zeros((nz,nz), dtype=img.dtype)

        img_proj = np.concatenate([
            np.concatenate([img_proj_xy, img_proj_yz.T], axis=1),
            np.concatenate([img_proj_xz,    img_offset], axis=1)], axis=0)

        return img_proj

    def get_stack_from_series_id(czi_path, channel, series_id, dim):

        # Extraction a series from multi-position CZI

        img_raw = []
        for slice in range(dim[0]):
            img_raw.append(bf.load_image(path=czi_path, series=series_id, c=channel, z=slice, t=0, rescale=False))
        img_raw = np.array(img_raw).reshape(*dim)

        return img_raw

    #
    # If feature extraction mode
    #

    if args["mode"] == "feature":

        # Specific imports

        from skimage.measure import label
        from skimage.transform import resize

        import javabridge as jv, bioformats as bf
        jv.start_vm(class_path=bf.JARS, max_heap_size='4G')

        from aicsfeature.extractor import cell, dna, structure

        #
        # Loading configuration form JSON file
        #

        with open(args["config"]+".json", "r") as fjson:
            config_full = json.load(fjson)
        config_name = config_full["name"]
        config_json = config_full["data"]
        print("Processing dataset:",config_name)

        # Checking whether static/imgs exist

        if not os.path.isdir("../engine/app/static/imgs/"):
            os.makedirs("../engine/app/static/imgs/")

        # Removing log file

        try:
            os.remove("log.json")
        except: pass

        # For each CZI

        print("Number of CZI files found:",len(config_json))
        
        df_meta = pd.DataFrame([])
        df_feat_dna = pd.DataFrame([])
        df_feat_mem = pd.DataFrame([])
        df_feat_str = pd.DataFrame([])

        for config_czi in config_json:

            czi_path = os.path.join(config_czi["raw_path"],config_czi["raw_name"])
            pxl_info = bf.OMEXML(bf.get_omexml_metadata(czi_path)).image(0).Pixels
            nx = pxl_info.SizeX
            ny = pxl_info.SizeY
            nz = pxl_info.SizeZ
            pixel_size_xy = config_czi["pixel_size_xy"]
            pixel_size_z = config_czi["pixel_size_z"]

            print("CZI:", config_czi["raw_name"], (nz,ny,nx))

            # For each series

            for position in config_czi["position"]:


                # Images from 880 microscope are single position, while
                # images from ZSD is multiposition and require us to
                # parse the seg file name to find the series id
                
                print("\tSeries:",position["name"])

                if config_czi["modality"] == 880:

                    series_id = 0

                else:

                    series_id = position["name"]
                    loc = series_id.find("Scene")

                    # Parse series ID in the segmentation name
                    # First positiuon is 0, alhtough it is shown as Scene-01
                    series_id = int(series_id.split("Scene")[1].split("-")[1]) - 1

                # Structure channel

                str_channel = 0
                if config_czi["minipipeline"] == "yes":
                    str_channel = 1
                    mem_channel = 0
                    dna_channel = 2

                img_raw_ok = False

                try:

                    img_str_raw = get_stack_from_series_id(czi_path=czi_path, channel=str_channel, series_id=series_id, dim=(nz,ny,nx))

                    img_mem_raw = get_stack_from_series_id(czi_path=czi_path, channel=mem_channel, series_id=series_id, dim=(nz,ny,nx))

                    img_dna_raw = get_stack_from_series_id(czi_path=czi_path, channel=dna_channel, series_id=series_id, dim=(nz,ny,nx))

                    img_raw_ok = True

                except:

                    # Log if not possible to read the raw data

                    with open("log.json", "a") as flog:
                        json.dump({
                            "czi_path": czi_path,
                            "position_name": position["name"],
                            "channel": str_channel,
                            "series_id": series_id,
                            "dim": (nz,ny,nz)}, flog, indent=4)
                    pass

                if img_raw_ok:

                    # Load segmentation & clear 1st and last slice

                    dna_seg_path = os.path.join(config_czi["seg_path"],position["name"])
                    mem_seg_path = os.path.join(config_czi["seg_path"],position["name"].replace("_nuc","_cell"))

                    img_dna_seg_all = skio.imread(dna_seg_path)
                    img_dna_seg_all[ 0,:,:] = 0
                    img_dna_seg_all[-1,:,:] = 0

                    img_mem_seg_all = skio.imread(mem_seg_path)
                    img_mem_seg_all[ 0,:,:] = 0
                    img_mem_seg_all[-1,:,:] = 0

                    # 880 images have been downsampled for segmentation

                    if config_czi["xy_upsample_seg"] > 1.0:

                        print("\tUpsampling segmentation, original:", img_dna_seg_all.shape, ", target:", nz,ny,ny, ", factor reported:", config_czi["xy_upsample_seg"])

                        img_dna_seg_all_dtype = img_dna_seg_all.dtype
                        img_dna_seg_all = resize(image = img_dna_seg_all,
                            output_shape = (nz, ny, nx),
                            order = 0,
                            preserve_range = True,
                            anti_aliasing = False,
                            mode = "constant")
                        img_dna_seg_all = img_dna_seg_all.astype(img_dna_seg_all_dtype)

                    # Analyze each cell

                    print("\tProcessing cell:")

                    for cell_id in position["cell_id"]:

                        print("\t\t",cell_id)

                        cell_name = position["name"].replace(".czi_nucWholeIndexSegScale.tiff","_cid_"+str(cell_id))

                        img_dna_seg = img_dna_seg_all.copy()
                        img_dna_seg[img_dna_seg!=cell_id] = 0
                        img_dna_seg[img_dna_seg==cell_id] = 1

                        img_mem_seg = img_mem_seg_all.copy()
                        img_mem_seg[img_mem_seg!=cell_id] = 0
                        img_mem_seg[img_mem_seg==cell_id] = 1

                        print("\t\t\tTesting object connectivity")

                        # Testing whether the nucleus has a unique component
                        # Excluding pixels=0 during bincount

                        img_dna_seg = label(img_dna_seg)

                        if img_dna_seg.max() > 1:
                            print("\t\t\t\tKeeping largest cc")
                            largest_cc = 1 + np.argmax(np.bincount(img_dna_seg.flat)[1:])
                            img_dna_seg[img_dna_seg!=largest_cc] = 0
                            img_dna_seg[img_dna_seg>0] = 1

                        img_mem_seg = label(img_mem_seg)

                        if img_mem_seg.max() > 1:
                            print("\t\t\t\tKeeping largest cc")
                            largest_cc = 1 + np.argmax(np.bincount(img_mem_seg.flat)[1:])
                            img_mem_seg[img_mem_seg!=largest_cc] = 0
                            img_mem_seg[img_mem_seg>0] = 1

                        # Cropping the nucleus

                        print("\t\t\tCropping ROI")

                        pxl_z, pxl_y, pxl_x = np.nonzero(img_mem_seg)

                        img_mem_seg_crop = img_mem_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                        img_mem_raw_crop = img_mem_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                        img_dna_seg_crop = img_dna_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                        img_dna_raw_crop = img_dna_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                        img_str_seg_crop = img_dna_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                        img_str_raw_crop = img_str_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                        img_dna_input = img_dna_seg_crop * img_dna_raw_crop

                        img_mem_input = img_mem_seg_crop * img_mem_raw_crop

                        img_str_input = img_str_seg_crop * img_str_raw_crop

                        # Rescale to isotropic volume

                        print("\t\t\tRescale z direction")

                        dim_z, dim_y, dim_x = img_dna_input.shape
                        dim_z = np.int((pixel_size_z/pixel_size_xy)*dim_z)

                        img_dna_input = resize(image=img_dna_input, output_shape=(dim_z,dim_y,dim_x), preserve_range=True, anti_aliasing=True, mode="constant")
                        img_dna_input = img_dna_input.astype(np.int64)

                        img_mem_input = resize(image=img_mem_input, output_shape=(dim_z,dim_y,dim_x), preserve_range=True, anti_aliasing=True, mode="constant")
                        img_mem_input = img_mem_input.astype(np.int64)

                        img_str_input = resize(image=img_str_input, output_shape=(dim_z,dim_y,dim_x), preserve_range=True, anti_aliasing=True, mode="constant")
                        img_str_input = img_str_input.astype(np.int64)

                        # Save the image for interactive view

                        print("\t\t\tSaving orthogonal projections:", cell_name)

                        img_png = get_orthogonal_projs(img_mem_input)
                        skio.imsave(os.path.join("../engine/app/static/imgs",cell_name+".png"),img_png.astype(np.uint16))

                        # Feature extraction

                        print("\t\t\tExtracting features...")

                        feat_dna = dna.get_features(img=img_dna_input, extra_features=["roundness"])
                        feat_dna["cell_id"] = cell_name

                        feat_mem = cell.get_features(img=img_mem_input)
                        feat_mem["cell_id"] = cell_name

                        feat_str = dna.get_features(img=img_str_input, extra_features=["io_intensity", "bright_spots", "roundness"])
                        feat_str["cell_id"] = cell_name

                        # Metadata

                        meta = pd.DataFrame({
                            "cell_id": cell_name,
                            "czi": config_czi["raw_name"],
                            "series_id": series_id,
                            "cell_seg_id": cell_id,
                            "condition": config_czi["condition"],
                            "cell_type": config_czi["cell_type"],
                            "structure": config_czi["structure"],
                            "scope": config_czi["modality"],
                            "minipipeline": config_czi["minipipeline"]}, index=[0])

                        df_meta = pd.concat([df_meta,meta], axis=0, ignore_index=True)
                        df_feat_mem = pd.concat([df_feat_mem,feat_mem], axis=0, ignore_index=True)
                        df_feat_dna = pd.concat([df_feat_dna,feat_dna], axis=0, ignore_index=True)
                        df_feat_str = pd.concat([df_feat_str,feat_str], axis=0, ignore_index=True)

                        df_meta.to_csv(os.path.join("../data-raw/",config_name+"_meta.csv"), index=False)
                        df_feat_mem.to_csv(os.path.join("../data-raw/",config_name+"_mem_feature.csv"), index=False)
                        df_feat_dna.to_csv(os.path.join("../data-raw/",config_name+"_dna_feature.csv"), index=False)
                        df_feat_str.to_csv(os.path.join("../data-raw/",config_name+"_str_feature.csv"), index=False)

        print("Done!")

    #
    # If config mode
    #

    if args["mode"] == "config":

        import glob

        # Removing existing config file

        try:
            os.remove("config.json")
        except: pass

        #
        # Load CSV file based on Jinaxu segmentation report in
        # /allen/aics/assay-dev/MicroscopyOtherData/Jianxu/Nucleus/nucleus_meta.csv
        #

        df = pd.read_csv(os.path.join("dataset/",args["dataset"]+".csv"), sep=";")

        position = []

        #
        # For each multiposition image
        #

        for index, czi in df.iterrows():

            print("Loading ", czi["raw_name"], "...")

            #
            # Checking if the CZI file exist
            #

            czi_path = os.path.join(czi["raw_path"], czi["raw_name"])

            if not os.path.isfile(czi_path):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.join(czi["raw_path"], czi["raw_name"]))

            #
            # Looking for corresponding segmentation files
            #

            seg_path = czi["seg_path"]

            filenames = glob.glob1(seg_path, czi["raw_name"].replace(".czi","") + "*nucWholeIndexImageScale.tiff")
            filenames.sort()

            #
            # For each position
            #

            list_of_ome_tif_files = []

            for fname in filenames:

                #
                # List cell ids here
                #

                img = skio.imread(os.path.join(seg_path,fname))

                cell_id = np.unique(img[img>0])

                list_of_ome_tif_files.append({"name": fname, "cell_id": cell_id.tolist()})

            position.append(list_of_ome_tif_files)

            if len(list_of_ome_tif_files) > 0:

                print("\tNumber of files found:", len(list_of_ome_tif_files))

            else:

                print("\tNo segmentation found.")

        df["position"] = position

        with open(args["dataset"]+".json", "w") as fj:
            json.dump({"name": args["dataset"], "data": df.to_dict("records")}, fj, indent=4)

    #
    # If contrast mode
    #

    if args["mode"] == "contrast":

        import glob

        #
        # Load processed CSV. This CSV does not include outliers.
        #

        df = pd.read_csv(os.path.join("../engine/data-processed/",args["dataset"]+".csv"), sep=";")

        #
        # For each cell
        #

        TARGET_MIN = 0.05;
        TARGET_MAX = 0.95;

        size_x, size_y = [], []
        percentile_min, percentile_max = [], []

        for index, cell in df.iterrows():

            img = skio.imread(os.path.join("../engine/app/static/imgs/",cell["cell_id"]+".png"))

            size_x.append(img.shape[1])
            size_y.append(img.shape[0])

            img = img[img>0].reshape(-1)

            percentile_min.append(np.percentile(img, q=100*TARGET_MIN))
            percentile_max.append(np.percentile(img, q=100*TARGET_MAX))

        pct_min = np.min(percentile_min)
        pct_max = np.max(percentile_max)

        print("Percentiles:", pct_min, pct_max)

        stack = np.zeros((len(size_y), np.max(size_y)+5, np.max(size_x)+5), dtype=np.uint8)

        #
        # Adjust contrast and save as jpg
        #

        for index, cell in df.iterrows():

            img = skio.imread(os.path.join("../engine/app/static/imgs/",cell["cell_id"]+".png"))

            img[(img>0)&(img<pct_min)] = pct_min
            img[(img>0)&(img>pct_max)] = pct_max

            img[img>0] = (255. * (img[img>0]-pct_min) / (pct_max-pct_min)).astype(np.uint8)
            img[0,0] = 255

            stack[index,:img.shape[0],:img.shape[1]] = img

            skio.imsave(os.path.join("../engine/app/static/imgs/",cell["cell_id"]+".jpg"), img)

        #skio.imsave(os.path.join("../engine/app/static/imgs/stack.tif"), stack)
