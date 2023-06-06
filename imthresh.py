import casatools
import glob

print("image thresholding: ")

# get list of model images
model_list = glob.glob("*.model.tt0")
ia = casatools.image()

for model in model_list:
    ttn_list = glob.glob(f"{model[:-1]}*")

    ia.open(model)
    im = ia.getchunk()
    ia.close()

    # make a mask of negative pixels
    mask = im < 1e-12

    # replace negative pixels with zeros in all taylor-term models
    for ttn in ttn_list:
        print(ttn)
        ia.open(ttn)
        imttn = ia.getchunk()
        imttn[mask] = 0.0
        ia.putchunk(imttn)
        ia.close()

