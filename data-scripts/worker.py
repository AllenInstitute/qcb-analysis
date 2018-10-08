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
                                            image (save cell images in static folder)", required=True)
    parser.add_argument("-c", "--config", help="Path to config json", required=True)
    args = vars(parser.parse_args())

    """
        Loading configuration form JSON file
    """

    with open(args["config"], "r") as fjson:
        config_json = json.load(fjson)

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

        from aicsfeature.extractor import mem, dna, structure

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

            #
            # Finding all cell with that structure
            #

            print("\nImage Type:[", struct["structure_name"], "]\n")

            if struct["structure_name"] not in ["mem", "dna"]:

                df_meta_struct = df_meta.loc[df_meta.structure_name==struct["structure_name"]]

            else:

                df_meta_struct = df_meta.copy()

            df_features = pd.DataFrame([])
            for row in tqdm(range(df_meta_struct.shape[0])):
                cid = df_meta_struct.index[row]
                seg_path = os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"])
                raw_path = os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"])
                SEG = skio.imread(seg_path)
                RAW = skio.imread(raw_path)

                #
                # Feature extraction for each cell
                #

                if struct["structure_name"] is "mem":

                    df_features = df_features.append(mem.GetFeatures(None,seg=RAW[1,:,:,:]*SEG[1,:,:,:]),
                                ignore_index=True)

                elif struct["structure_name"] is "dna":

                    df_features = df_features.append(dna.GetFeatures(None,seg=RAW[0,:,:,:]*SEG[0,:,:,:]),
                                ignore_index=True)
                else:

                    df_features = df_features.append(structure.GetFeatures(None,seg=RAW[2,:,:,:]*SEG[2,:,:,:],
                                extra_features=struct["extra_features"]),
                                ignore_index=True)

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

        """
            Creates an image with custom colors
        """

        def fix_cell_color(img):
            for i in [2,0,1]:
                for j in [2,0,1]:
                    if i != j:
                        img[img[:,:,i]>0,j] = 0
            img[np.all(img==0,axis=2),:] = 255
            return img

        """
            Load a given cell
        """

        def get_cell_image(cell_img_path):
            img = io.imread(cell_img_path)
            img = np.max(img,axis=1)
            img = np.swapaxes(img,0,2)
            img = img[:,:,:3]
            img = fixColor(img)
            return img

        for row in tqdm(range(df_meta.shape[0])):
            cid = df_meta_struct.index[row]
            img = get_cell_image(os.path.join(config_json["cell_info"],cid,config_json["seg_prefix"])
            skio.imsave(os.path.join('../engine/app/static/imgs',cid+'.jpg'),img)
