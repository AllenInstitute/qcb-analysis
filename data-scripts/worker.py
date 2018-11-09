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
    parser.add_argument("-m", "--mode", help="download (download metadata), \
                                            feature (feature extraction) \
                                            image (save cell images in static folder)\
                                            check (check dataframe produced by feature extraction)\
                                            process (process tables for further analyzes)", required=True)
    parser.add_argument("-c", "--config", help="Path to config json", required=True)
    parser.add_argument("-s", "--start", nargs="?", type=int, default=0, const=0, required=False)
    args = vars(parser.parse_args())

    #
    # Loading configuration form JSON file
    #

    with open(args["config"], "r") as fjson:
        config_json = json.load(fjson)

    # mem_ch = config_json["mem_channel"]
    # dna_ch = config_json["dna_channel"]
    # str_ch = config_json["str_channel"]

    #
    # If download mode
    #

    if args["mode"] == "download":

        import datasetdatabase as dsdb

        #
        # Multiprocessing
        #

        os.environ["DSDB_PROCESS_LIMIT"] = "16"

        #
        # Loading metadata from database
        #

        prod = dsdb.DatasetDatabase(config=config_json["database"])
        ds_meta = prod.get_dataset(name=config_json["meta"])

        with open(os.path.join("../data-raw/",config_json["meta"]+".pkl"), "wb") as fp:
            pickle.dump(ds_meta.ds,fp)

    #
    # If feature extraction mode
    #

    if args["mode"] == "feature":

        # Specific imports

        import javabridge as jv, bioformats as bf
        jv.start_vm(class_path=bf.JARS, max_heap_size='4G')

        from aicsfeature.extractor import cell, dna, structure

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

        # For each CZI

        print("Number of CZI files found:",len(config_json))
        
        Table = pd.DataFrame([])

        for config_czi in config_json:

            czi_path = os.path.join(config_czi["raw_path"],config_czi["raw_name"])
            pxl_info = bf.OMEXML(bf.get_omexml_metadata(czi_path)).image(0).Pixels
            nx = pxl_info.SizeX
            ny = pxl_info.SizeY
            nz = pxl_info.SizeZ

            print("CZI:", config_czi["raw_name"], (nz,ny,nx))

            # For each series

            for position in config_czi["ids_cell"]:

                # Parse series ID in the segmentation name
                # First positiuon is 0, alhtough it is shown as Scene-01

                print("\tSeries:",position["name"])

                series_id = position["name"]
                loc = series_id.find("Scene")

                # First positiuon is 0, alhtough it is shown as Scene-01
                series_id = np.int(series_id[(loc+6):(loc+8)]) - 1

                img_raw = get_stack_from_series_id(czi_path=czi_path, channel=config_czi["dna_channel"], series_id=series_id, dim=(nz,ny,nx))

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

                    pxl_z, pxl_y, pxl_x = np.nonzero(img_seg)

                    img_seg_crop = img_seg[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]
                    img_raw_crop = img_raw[pxl_z.min():(pxl_z.max()+1),pxl_y.min():(pxl_y.max()+1),pxl_x.min():(pxl_x.max()+1)]

                    img_input = img_seg_crop*img_raw_crop

                    # Save the image for interactive view

                    skio.imsave(os.path.join('../engine/app/static/imgs',cell_name+'.jpg'),img_input.max(axis=0))

                    # Feature extraction

                    df = dna.get_features(img=img_input, extra_features=["io_intensity", "bright_spots"])

                    # >>>>>>> split df in meta and features

                    # Metadata

                    df["czi"] = config_czi["raw_name"]
                    df["series_id"] = series_id
                    df["cell_id"] = cell_id
                    df["condition"] = config_czi["czi_label"]

                    print("\t\tSeries:", series_id, "Cell:", cell_id)

                    Table = pd.concat([Table,df], axis=0, ignore_index=True)

                    Table.to_csv("result.csv", index=False)

        print("Done!")


    #
    # If image mode
    #

    if args["mode"] == "image":

        #
        # Loading metadata
        #

        with open(os.path.join("../data-raw/",config_json["meta"]+".pkl"), "rb") as fp:
            df_meta = pickle.load(fp)

        df_meta = df_meta.set_index(config_json["id"])

        #
        # Checking whether static/imgs exist
        #

        if not os.path.isdir("../engine/app/static/imgs/"):
            os.makedirs("../engine/app/static/imgs/")

        #
        # Creates an image with custom colors
        #

        def fix_cell_color(img):
            for i in [str_ch,dna_ch,mem_ch]:
                for j in [str_ch,dna_ch,mem_ch]:
                    if i != j:
                        img[img[:,:,i]>0,j] = 0
            img[np.all(img==0,axis=str_ch),:] = 255
            return img

        def get_cell_image(cell_img_path):
            img = skio.imread(cell_img_path)
            img = np.max(img,axis=mem_ch)
            img = np.swapaxes(img,dna_ch,str_ch)
            img = img[:,:,:3]
            img = fix_cell_color(img)
            return img

        for row in tqdm(range(df_meta.shape[0])):
            cid = df_meta.index[row]
            img = get_cell_image(os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"]))
            skio.imsave(os.path.join('../engine/app/static/imgs',cid+'.jpg'),img)

    #
    # If check mode
    #

    if args["mode"] == "check":

        import time

        #
        # For each structure in the config file
        #

        report = pd.DataFrame([])
        for struct in config_json["run"]:

            #
            # Finding all cell with that structure
            #

            mod_date = time.ctime(os.path.getmtime(os.path.join("../data-raw/",struct["save_as"])))

            with open(os.path.join("../data-raw/",struct["save_as"]), "rb") as fp:
                df = pickle.load(fp)

            print("\n## "+struct["save_as"]+" ##\n")
            print(df.head())

            report = report.append({"name": struct["save_as"].replace(".pkl",""), "rows": df.shape[0], "cols": df.shape[1], "modified": mod_date}, sort=True, ignore_index=True)

        report[["cols","rows"]] = report[["cols","rows"]].astype(np.int)
        print(report[["name","cols","rows","modified"]])

    #
    # If process mode
    #

    if args["mode"] == "process":

        #
        # Loading metadata
        #

        with open(os.path.join("../data-raw/",config_json["meta"]+".pkl"), "rb") as fp:
            df_full = pickle.load(fp)

        df_full = df_full.set_index(config_json["id"])

        #
        # Loading features
        #

        for struct in config_json["run"]:

            if struct["status"] == "on":

                #
                # Finding all cell with that structure
                #

                with open(os.path.join("../data-raw/",struct["save_as"]), "rb") as fp:
                    df_str_fea = pickle.load(fp)
                    df_full = df_full.join(df_str_fea)

        df_full.to_csv("../engine/data-processed/data.csv")
