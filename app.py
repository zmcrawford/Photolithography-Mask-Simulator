import os
from datetime import datetime
from flask import Flask, flash, redirect, render_template, request, jsonify
import numpy as np
import base64
from io import BytesIO
from PIL import Image
from cffi import FFI

# https://cffi.readthedocs.io/en/stable/overview.html#main-mode-of-usage
ffi = FFI()
ffi.cdef("""
void compute_pixel(const uint8_t* input,
                   uint8_t* output,
                   const int width,
                   const int height,
                   float k1,
                   float NA,
                   float wavelength,
                   int psf_type);
""")

# Load the precompiled shared library
app = Flask(__name__)

lib = ffi.dlopen("./static/psf_sim.so")


@app.route("/")
def index():

    return render_template("index.html")


@app.route("/simulate", methods=['POST'])
def simulate():
    data = request.get_json()
    wavelength = float(data['wavelength'])
    NA = float(data['NA'])
    k1 = float(data['k1'])
    if data['method'] == 'gauss':
        psf_method = 0
    elif data['method'] == 'airy':
        psf_method = 1
    # TODO implement apology/else/safety things

    width = 500
    height = 500
    img_base64 = data["image"].split(",")[1]  # remove prefix in json

    # conversions adapted from (https://jdhao.github.io/2020/03/17/base64_opencv_pil_image_conversion)
    img_bytes = base64.b64decode(img_base64)

    # Open as PIL Image
    image = Image.open(BytesIO(img_bytes))

    # convert image to 1d binary image
    image = image.convert("1")
    binary_mask_1d = np.array(image, dtype=np.uint8).flatten()  # 1d array, 1 byte per pixel
    output_image = np.zeros_like(binary_mask_1d, dtype=np.uint8)

    input_c = ffi.cast("uint8_t *",  binary_mask_1d.ctypes.data)
    output_c = ffi.cast("uint8_t *", output_image.ctypes.data)

    # pass to precompiled c library
    lib.compute_pixel(input_c, output_c, width, height, k1, NA, wavelength, psf_method)

    # Wrap output correctly
    output_np = np.frombuffer(ffi.buffer(output_c, height*width), dtype=np.uint8)
    output_np = output_np.reshape((height, width))

    # Convert to PIL and save
    img = Image.fromarray(output_np)
    img.save('static/output.png')
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return jsonify({"image": f"data:image/png;base64,{img_str}"})
