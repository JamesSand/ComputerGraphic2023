
from PIL import Image
import argparse
ratio = 0.5
parser = argparse.ArgumentParser()
parser.add_argument("--r", type=str)
parser.add_argument("--g", type=str)
parser.add_argument("--b", type=str)
parser.add_argument("--out", type=str)
args = parser.parse_args()
aImg, bImg, cImg = Image.open(args.r), Image.open(args.g), Image.open(args.b)
bImg=Image.blend(aImg,bImg,ratio)
out=Image.blend(bImg,cImg,ratio)
out.save(args.out)

