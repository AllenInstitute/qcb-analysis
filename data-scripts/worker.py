import os
import json
import pickle
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import io as skio

if __name__ == "__main__":

    """
        External parameters
    """

    parser = argparse.ArgumentParser(description="Interface for downloading/exploring data from database and extracting features for QCB")
    parser.add_argument("-m", "--mode", help="download (download metadata), \
                                            feature (feature extraction) \
                                            image (save cell images in static folder)\
                                            check (check dataframe produced by feature extraction)\
                                            process (process tables for further analyzes)", required=True)
    parser.add_argument("-c", "--config", help="Path to config json", required=True)
    parser.add_argument("-s", "--start", nargs="?", type=int, default=0, const=0, required=False)
    args = vars(parser.parse_args())

    """
        Loading configuration form JSON file
    """

    with open(args["config"], "r") as fjson:
        config_json = json.load(fjson)

    mem_ch = config_json["mem_channel"]
    dna_ch = config_json["dna_channel"]
    str_ch = config_json["str_channel"]

    """
        If download mode
    """

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

    """
        If feature extraction mode
    """

    if args["mode"] == "feature":

        from aicsfeature.extractor import cell, dna, structure

        #
        # Loading metadata
        #

        with open(os.path.join("../data-raw/",config_json["meta"]+".pkl"), "rb") as fp:
            df_meta = pickle.load(fp)

        #
        # Indexing metadata dataframe by id specified in config_json. Keep in mind
        # that this id will not be maintained in database and a explicity id column
        # is required to be created before uploading a dataset.
        #

        df_meta = df_meta.set_index(config_json["id"])

        #
        # For each structure in the config file
        #

        for struct in config_json["run"]:


            if struct["status"] == "on":
            
                #
                # Finding all cell with that structure
                #

                print("\nImage Type:[", struct["structure_name"], "]\n")

                if struct["structure_name"] not in ["cell", "dna"]:

                    df_meta_struct = df_meta.loc[df_meta.structure_name==struct["structure_name"]]

                else:

                    df_meta_struct = df_meta.copy()

                df_features = pd.DataFrame([])
                for row in tqdm(range(df_meta_struct.shape[0])):

                    if row >= args["start"]:

                        cid = df_meta_struct.index[row]
                        seg_path = os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"])
                        raw_path = os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"])
                        SEG = skio.imread(seg_path)
                        RAW = skio.imread(raw_path)

                        #
                        # Feature extraction for each cell
                        #

                        if struct["structure_name"] == "cell":
                            df_features = df_features.append(cell.get_features(img=RAW[mem_ch,:,:,:]*SEG[mem_ch,:,:,:]),
                                ignore_index=True, sort=True)

                        elif struct["structure_name"] == "dna":
                            df_features = df_features.append(dna.get_features(img=RAW[dna_ch,:,:,:]*SEG[dna_ch,:,:,:]),
                                ignore_index=True, sort=True)
                        else:
                            df_features = df_features.append(structure.get_features(img=RAW[str_ch,:,:,:]*SEG[str_ch,:,:,:],
                                extra_features=struct["extra_features"]),
                                ignore_index=True, sort=True)
            
                df_features.index = df_meta_struct.index

                #
                # Save as pickle
                #

                with open(os.path.join("../data-raw/",struct["save_as"]), "wb") as fp:
                    pickle.dump(df_features,fp)

    """
        If image mode
    """

    if args["mode"] == "image":

        #
        # Loading metadata
        #

        with open(os.path.join("../data-raw/",config_json["meta"]+".pkl"), "rb") as fp:
            df_meta = pickle.load(fp)

        df_meta = df_meta.set_index(config_json["id"])

        """
            Checking whether static/imgs exist
        """

        if not os.path.isdir("../engine/app/static/imgs/"):
            os.makedirs("../engine/app/static/imgs/")

        """
            Creates an image with custom colors
        """

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

    """
        If check mode
    """

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

    """
        If process mode
    """

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
