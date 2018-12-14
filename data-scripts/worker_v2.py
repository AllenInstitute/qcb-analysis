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

from worker_tools import make_pdf_report

if __name__ == "__main__":

    # External parameters

    parser = argparse.ArgumentParser(description="Interface for downloading/exploring data from database and extracting features for QCB")
    parser.add_argument("-m", "--mode", help="feature (feature extraction) \
                                            contrast (adjust the contrast of images of a particular dataset)\
                                            config (generate the config json file based on a csv table)", required=True)
    parser.add_argument("-c", "--config", nargs="?", type=str, default="", const="", required=False)
    parser.add_argument("-d", "--dataset", nargs="?", type=str, default="", const="", required=False)
    args = vars(parser.parse_args())

    # Auxiliar functions

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

    def get_image_dimension(czi_path):

        '''
            This function assumes that all images in the CZI file
            have same dimension and therefore, we only need to
            retrieve the dimensions of the 1st image.
        '''

        pxl_info = bf.OMEXML(bf.get_omexml_metadata(czi_path)).image(0).Pixels
        nx = pxl_info.SizeX
        ny = pxl_info.SizeY
        nz = pxl_info.SizeZ

        return (nz, ny, nx)

    def get_stack_from_series_id(czi_path, channel, series_id, dim):

        '''
            Extracts a z-stack (position) from multi-position CZI
            Slice-by-slice is required in the current bio-format
            package
        '''

        img_raw = []
        for slice in range(dim[0]):
            img_raw.append(bf.load_image(path=czi_path, series=series_id, c=channel, z=slice, t=0, rescale=False))
        img_raw = np.array(img_raw).reshape(*dim)

        return img_raw

    def get_list_valid_positions(text):

        '''
            Valid positions can be provided as list [1,2,3...] or
            sequence 1:10, or a single position. Or it can also be
            a string:
                force0 - it is always going to map to first position
                any - any position is valid
        '''

        list_of_positions = None

        if str(text).isdigit():

            list_of_positions = [int(text)]
        
        elif ":" in text:

            i = np.int(text.split(":")[0])
            j = np.int(text.split(":")[1])

            list_of_positions = np.arange(i,j+1,1).tolist()

        elif "," in text:

            list_of_positions = list(map(int, text.split(",")))

        elif "force0" == text:

            list_of_positions = -1

        return list_of_positions

    def print_log(text, log_file):
        print(text)
        with open(log_file, "a") as flog:
            flog.write(text+"\n")

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
