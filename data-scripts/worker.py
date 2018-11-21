import os
import json
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
    parser.add_argument("-d", "--dataset", nargs="?", type=str, default="", const="", required=False)
    parser.add_argument("-s", "--start", nargs="?", type=int, default=0, const=0, required=False)
    args = vars(parser.parse_args())

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

        with open("config.json", "r") as fjson:
            config_full = json.load(fjson)
        config_name = config_full["name"]
        config_json = config_full["data"]
        print("Processing dataset:",config_name)

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
        df_feat = pd.DataFrame([])

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

                # Parse series ID in the segmentation name
                # First positiuon is 0, alhtough it is shown as Scene-01

                print("\tSeries:",position["name"])

                series_id = position["name"]
                loc = series_id.find("Scene")

                # First positiuon is 0, alhtough it is shown as Scene-01
                series_id = int(series_id.split("Scene")[1].split("-")[1]) - 1

                # Structure channel

                str_channel = 0
                if config_czi["minipipeline"] == "yes":
                    str_channel = 1

                img_raw_ok = False

                try:

                    img_raw = get_stack_from_series_id(czi_path=czi_path, channel=str_channel, series_id=series_id, dim=(nz,ny,nx))

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

                    seg_path = os.path.join(config_czi["seg_path"],position["name"])
                    img_seg_all = skio.imread(seg_path)
                    img_seg_all[ 0,:,:] = 0
                    img_seg_all[-1,:,:] = 0

                    # Analyze each cell

                    for cell_id in position["cell_id"]:

                        cell_name = position["name"].replace(".ome.tif","_cid_"+str(cell_id))

                        img_seg = img_seg_all.copy()
                        img_seg[img_seg!=cell_id] = 0
                        img_seg[img_seg==cell_id] = 1

                        # Testing whether the nucleus has a unique component
                        # Excluding pixels=0 during bincount

                        img_seg = label(img_seg)

                        if img_seg.max() > 1:
                            largest_cc = 1 + np.argmax(np.bincount(img_seg.flat)[1:])
                            img_seg[img_seg!=largest_cc] = 0
                            img_seg[img_seg>0] = 1

                        # Cropping the nucleus

                        pxl_z, pxl_y, pxl_x = np.nonzero(img_seg)

                        img_seg_crop = img_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                        img_raw_crop = img_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                        img_input = img_seg_crop * img_raw_crop

                        # Rescale to isotropic volume

                        dim_z, dim_y, dim_x = img_input.shape
                        dim_z = np.int((pixel_size_z/pixel_size_xy)*dim_z)
                        img_input = resize(image=img_input, output_shape=(dim_z,dim_y,dim_x), preserve_range=True, anti_aliasing=True, mode="constant")
                        img_input = img_input.astype(np.int64)

                        # Save the image for interactive view

                        img_png = get_orthogonal_projs(img_input)
                        skio.imsave(os.path.join("../engine/app/static/imgs",cell_name+".png"),img_png.astype(np.uint16))

                        # Feature extraction

                        feat = dna.get_features(img=img_input, extra_features=["io_intensity", "bright_spots", "roundness"])
                        feat["cell_id"] = cell_name
                        
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

                        print("\t\tSeries:", series_id, "Cell:", cell_id)

                        df_meta = pd.concat([df_meta,meta], axis=0, ignore_index=True)
                        df_feat = pd.concat([df_feat,feat], axis=0, ignore_index=True)

                        df_meta.to_csv(os.path.join("../data-raw/",config_name+"_meta.csv"), index=False)
                        df_feat.to_csv(os.path.join("../data-raw/",config_name+"_feature.csv"), index=False)

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

            print("Loading ", czi["experiment_id"], "...")

            seg_path = czi["seg_path"]

            filenames = glob.glob1(seg_path,"*.ome.tif")
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

        df["position"] = position

        with open("config.json", "w") as fj:
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

        skio.imsave(os.path.join("../engine/app/static/imgs/stack.tif"), stack)