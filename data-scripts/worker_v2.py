import os
import json
import errno
import pickle
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio

from colorama import init, Fore
init(autoreset=True)

from worker_tools import *

if __name__ == "__main__":

    # External parameters

    parser = argparse.ArgumentParser(description="Interface for downloading/exploring data from database and extracting features for QCB")
    parser.add_argument("-m", "--mode", help="feature (feature extraction) \
                                            contrast (adjust the contrast of images of a particular dataset)\
                                            config (generate the config json file based on a csv table)", required=True)
    parser.add_argument("-c", "--config", nargs="?", type=str, default="", const="", required=False)
    parser.add_argument("-d", "--dataset", nargs="?", type=str, default="", const="", required=False)
    args = vars(parser.parse_args())

    #
    # If config mode
    #

    if args["mode"] == "config":

        import re
        import glob
        from skimage.transform import resize

        import javabridge as jv, bioformats as bf
        jv.start_vm(class_path=bf.JARS, max_heap_size='4G')

        # Load dataset CSV

        df = pd.read_csv(os.path.join("dataset/",args["dataset"]+".csv"), sep=";")

        # Logging

        log_file = args["dataset"]+".log"

        try:
            os.remove(log_file)
        except: pass

        # For each multiposition image

        position = []

        for index, czi in df.iterrows():

            # Checking if the CZI file exist

            czi_path = os.path.join(czi["czi_path"], czi["czi_name"]+".czi")

            print_log(":: "+czi["czi_name"], log_file=log_file)
            print_log(":: "+czi_path, log_file=log_file)

            if not os.path.isfile(czi_path):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), czi_path)

            # Calculating image dimension

            czi_dim = get_image_dimension(czi_path=czi_path)

            # Looking for corresponding segmentation files

            print_log(":: "+czi["seg_path"], log_file=log_file)

            filenames = glob.glob1(czi["seg_path"], czi["czi_name"] + "*" + czi["seg_suffix"] + "*")
            filenames.sort()

            valid_positions = get_list_valid_positions(czi["position"])

            # For each segmentation file

            list_of_seg_files = []

            for fname in filenames:

                #
                # List cell ids here
                #

                is_valid = True

                seg_path = os.path.join(czi["seg_path"], fname)

                position_id = int(re.search(r'Scene-\s*(\d+)', fname)[1])

                # List of valid positions

                if valid_positions is not None:

                    if valid_positions == -1:

                        position_id = 1

                    elif position_id not in valid_positions:

                        is_valid = False

                if is_valid:

                    # Open raw image

                    img_raw = get_stack_from_series_id(czi_path = czi_path,
                        channel = int(czi["channel"]),
                        series_id = position_id-1,
                        dim = czi_dim)

                    # Open corresponding segmentation

                    img_seg = skio.imread(seg_path)

                    cell_id = np.unique(img_seg[img_seg>0])

                    # Resize if necessary

                    if img_seg.shape != img_raw.shape:

                        print_log("\tUpsampling seg image", log_file=log_file)

                        img_seg_dtype = img_seg.dtype
                        img_seg = resize(image = img_seg,
                            output_shape = czi_dim,
                            order = 0,
                            preserve_range = True,
                            anti_aliasing = False,
                            mode = "constant")
                        img_seg = img_seg.astype(img_seg_dtype)

                    # Calculate background and foreground histograms

                    hist_foreground_y = np.bincount(img_raw[img_seg>=1])
                    hist_foreground_x = np.arange(0,hist_foreground_y.shape[0])
                    hist_foreground_x = hist_foreground_x[hist_foreground_y>1].tolist()
                    hist_foreground_y = hist_foreground_y[hist_foreground_y>1].tolist()

                    hist_background_y = np.bincount(img_raw[img_seg==0])
                    hist_background_x = np.arange(0,hist_background_y.shape[0])
                    hist_background_x = hist_background_x[hist_background_y>1].tolist()
                    hist_background_y = hist_background_y[hist_background_y>1].tolist()

                    list_of_seg_files.append({"name": fname,
                        "cell_id": cell_id.tolist(),
                        "hist_foreground": {"x": hist_foreground_x, "y": hist_foreground_y},
                        "hist_background": {"x": hist_background_x, "y": hist_background_y}})

                    print_log(Fore.GREEN+"\t[\u2713] "+fname, log_file=log_file)

                else:

                    print_log(Fore.RED+"\t[x] "+fname, log_file=log_file)

            position.append(list_of_seg_files)

            if len(list_of_seg_files) > 0:

                print_log("\t"+"Number of files found: "+str(len(list_of_seg_files)), log_file=log_file)

            else:

                print_log("\t"+"No segmentation found.", log_file=log_file)

        df["position"] = position

        # Save in JSON

        with open(args["dataset"]+".json", "w") as fj:
            json.dump({"name": args["dataset"], "data": df.to_dict("records")}, fj, indent=4)

    jv.kill_vm()

    # Create a report

    make_pdf_report(file_name=args["dataset"])

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


                # Images from 880 microscope are single position, while
                # images from ZSD is multiposition and require us to
                # parse the seg file name to find the series id
                
                print("\tSeries:",position["name"])

                if config_czi["modality"] == 880:

                    series_id = 0

                else:

                    series_id = position["name"]

                    # Parse series ID in the segmentation name
                    # First positiuon is 0, alhtough it is shown as Scene-01
                    if config_czi["force_single_scene"] == "yes":
                        series_id = 0
                    else:
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

                    # 880 images have been downsampled for segmentation

                    if config_czi["xy_upsample_seg"] > 1.0:

                        print("\tUpsampling segmentation, original:", img_seg_all.shape, ", target:", nz,ny,ny, ", factor reported:", config_czi["xy_upsample_seg"])

                        img_seg_all_dtype = img_seg_all.dtype
                        img_seg_all = resize(image = img_seg_all,
                            output_shape = (nz, ny, nx),
                            order = 0,
                            preserve_range = True,
                            anti_aliasing = False,
                            mode = "constant")
                        img_seg_all = img_seg_all.astype(img_seg_all_dtype)

                    # Analyze each cell

                    print("\tProcessing cell:")

                    for cell_id in position["cell_id"]:

                        print("\t\t",cell_id)

                        cell_name = position["name"].replace(".ome.tif","_cid_"+str(cell_id))

                        img_seg = img_seg_all.copy()
                        img_seg[img_seg!=cell_id] = 0
                        img_seg[img_seg==cell_id] = 1

                        print("\t\t\tTesting object connectivity")

                        # Testing whether the nucleus has a unique component
                        # Excluding pixels=0 during bincount

                        img_seg = label(img_seg)

                        if img_seg.max() > 1:
                            largest_cc = 1 + np.argmax(np.bincount(img_seg.flat)[1:])
                            img_seg[img_seg!=largest_cc] = 0
                            img_seg[img_seg>0] = 1

                        # Cropping the nucleus

                        print("\t\t\tCropping ROI")

                        pxl_z, pxl_y, pxl_x = np.nonzero(img_seg)

                        img_seg_crop = img_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                        img_raw_crop = img_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                        img_input = img_seg_crop * img_raw_crop

                        # Rescale to isotropic volume

                        print("\t\t\tRescale z direction")

                        dim_z, dim_y, dim_x = img_input.shape
                        dim_z = np.int((pixel_size_z/pixel_size_xy)*dim_z)
                        img_input = resize(image=img_input, output_shape=(dim_z,dim_y,dim_x), preserve_range=True, anti_aliasing=True, mode="constant")
                        img_input = img_input.astype(np.int64)

                        # Save the image for interactive view

                        print("\t\t\tSaving orthogonal projections:", cell_name)

                        img_png = get_orthogonal_projs(img_input)
                        skio.imsave(os.path.join("../engine/app/static/imgs",cell_name+".png"),img_png.astype(np.uint16))

                        # Feature extraction

                        print("\t\t\tExtracting features...")

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

                        df_meta = pd.concat([df_meta,meta], axis=0, ignore_index=True)
                        df_feat = pd.concat([df_feat,feat], axis=0, ignore_index=True)

                        df_meta.to_csv(os.path.join("../data-raw/",config_name+"_meta.csv"), index=False)
                        df_feat.to_csv(os.path.join("../data-raw/",config_name+"_feature.csv"), index=False)

        print("Done!")