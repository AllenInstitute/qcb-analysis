import json
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages

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

def make_pdf_report(file_name):

    with open(file_name+".json", "r") as fjson:
        config = json.load(fjson)

    plots = []
    for czi in config["data"]:

        fig, ax = plt.subplots(2,2,
                               figsize=(12, 6+6*len(czi["position"])/50.),
                               gridspec_kw={'height_ratios':[10, int(1+9*len(czi["position"])/50.)]})
        fig.suptitle(czi["czi_name"])

        inset = inset_axes(ax[0,0], width="40%", height="40%")
        
        names = []
        xmin_f, xmin_b = [], []
        xmax_f, xmax_b = [], []
        for pos_id, position in enumerate(czi["position"]):
            hf = position["hist_foreground"]
            hb = position["hist_background"]

            fx = np.array(hf["x"])
            fy = np.array(hf["y"])
            fy = fy / fy.sum()
            avg_f = (fx*fy).sum()
            std_f = np.sqrt(((fx**2)*fy).sum()-(avg_f**2))

            bx = np.array(hb["x"])
            by = np.array(hb["y"])
            by = by / by.sum()
            avg_b = (bx*by).sum()
            std_b = np.sqrt(((bx**2)*by).sum()-(avg_b**2))
            
            xmin_b.append(avg_b-2*std_b)
            xmax_b.append(avg_b+2*std_b)
            xmin_f.append(avg_f-2*std_f)
            xmax_f.append(avg_f+2*std_f)

            rcolor = np.random.rand(3,)
            ax[0,0].plot(bx, by, c=rcolor)
            inset.set_xticks([])
            inset.plot(pos_id, avg_b, "o", c=rcolor, markersize=4)
            inset.plot(pos_id, bx.min(), "x", c="black", markersize=2)
            
            ax[0,1].plot(fx, fy, c=rcolor)

            ax[1,1].barh(pos_id, len(position["cell_id"]), color=rcolor)
            
            names.append(position["name"])
            
        ax[0,0].set_xlim([np.min(xmin_b),np.max(xmax_b)])
        ax[0,1].set_xlim([np.min(xmin_f),np.max(xmax_f)])
        ax[1,0].axis("off")
        ax[1,1].yaxis.set_ticks(np.arange(0,len(names)))
        ax[1,1].set_yticklabels(names)

        plots.append(fig)

    pp = PdfPages(file_name+".pdf")
    for fig in plots:
        pp.savefig(fig, dpi = 300, transparent = True)
    pp.close()
