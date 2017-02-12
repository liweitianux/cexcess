#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Created: 2016-05-08
#

"""
Make a montage of the input 4 images by utilizing the `montage`
tool from `ImageMagick`.

Montage:
+----------+----------+
|          |          |
|  image1  |  image2  |
|          |          |
+----------+----------+
|          |          |
|  image3  |  image4  |
|          |          |
+----------+----------+

Note:
The resolution/dimension of the `image1` is taken as the reference, and
the other 3 images are re-sized to match it in order to make an compact
and grid-aligned montage.


Credits:
[1] ImageMagick v6 Examples -- Montage, Arrays of Images
    https://www.imagemagick.org/Usage/montage/
[2] ImageMagick montage: Resizing only one input image
    https://stackoverflow.com/a/26273173/4856091
"""


import subprocess
import argparse


def get_image_size(image):
    """
    Get the size (dimension) of the image by utilizing `identify`.
    """
    results = subprocess.run(["identify", "-ping", "-format", "%w %h", image],
                             check=True, stdout=subprocess.PIPE)
    width, height = map(int, results.stdout.decode("utf-8").split())
    return (width, height)


def montage_images(infiles, outfile):
    """
    Combine the 4 input images to make a montage by utilizing `montage`.
    """
    img1, img2, img3, img4 = infiles
    # `img1` is taken as the reference image (top-left)
    w_ref, h_ref = get_image_size(img1)
    w_img2, h_img2 = get_image_size(img2)
    w_img3, h_img3 = get_image_size(img3)
    # `img2` is re-sized to have the same height as `img1` (top-right)
    ratio_img2 = h_ref / h_img2
    W_img2, H_img2 = int(w_img2*ratio_img2), int(h_img2*ratio_img2)
    # `img3` is re-sized to match the width of `img1` (bottom-left)
    ratio_img3 = w_ref / w_img3
    W_img3, H_img3 = int(w_img3*ratio_img3), int(h_img3*ratio_img3)
    # `img4` is simply re-sized to fill the bottom-right area
    W_img4, H_img4 = W_img2, H_img3
    subprocess.run(args=[
        "montage", img1,
        "(", img2, "-resize", "%dx%d" % (W_img2, H_img2), ")",
        "(", img3, "-resize", "%dx%d" % (W_img3, H_img3), ")",
        "(", img4, "-resize", "%dx%d" % (W_img4, H_img4), ")",
        "-geometry", "+1+1", "-tile", "2x2", outfile,
    ])


def main():
    parser = argparse.ArgumentParser(
        description="Make a montage of the input 4 images")
    parser.add_argument("-o", "--outfile", dest="outfile", required=True,
                        help="output filename of the montage image")
    parser.add_argument("images", nargs=4,
                        help="filenames of 4 input images")
    args = parser.parse_args()

    montage_images(infiles=args.images, outfile=args.outfile)


if __name__ == "__main__":
    main()
