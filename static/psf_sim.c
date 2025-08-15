#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define GAUSSIAN 0
#define AIRY 1
#define CANVAS_WIDTH 500
#define MASK_WIDTH 1000 // 1000 nm = 1 um
#define PI 3.14159265358979323846

// Compute PSF pixel values
void compute_pixel(const uint8_t *input, uint8_t *output, const int width, const int height,
                   float k1, float NA, float wavelength, int psf_type)
{

    float pixel_size = (float) MASK_WIDTH / (float) CANVAS_WIDTH; // nm per pixel
    float rayleigh = k1 * wavelength / (NA);                      // Rayleigh criterion (nm)

    // apply gaussian psf
    if (psf_type == GAUSSIAN)
    {
        float sigma_px = rayleigh / pixel_size; // convert nm to pixels
        float *temp = malloc(width * height * sizeof(float));
        float *temp2 = malloc(width * height * sizeof(float));
        // build gaussian kernel (definition: G(x,y) = 1/(norm constant)exp(-(x^2+y^2)/(2\sigma^2)),
        // (https://staff.fnwi.uva.nl/r.vandenboomgaard/ComputerVision/LectureNotes/IP/LocalStructure/GaussianDerivatives.html)
        int kernel_size = (int) (sigma_px * 6); // +- 3 sigma
        if (kernel_size % 2 == 0)
            kernel_size++; // need an odd kernel size
        float *kernel = malloc(kernel_size * kernel_size * sizeof(float));
        float sum = 0.0f;
        int half = kernel_size / 2; // distance from center of kernel

        // make 1D kernel which I'll apply horizontally then vertically to avoid 4x nested for loops
        // (matrix is outer prod of 2 1d vectors, nicely decomposes due to symmetry)
        for (int i = 0; i < kernel_size; i++)
        {
            kernel[i] = exp(-((i - half) * (i - half)) / (2 * sigma_px * sigma_px));
            sum += kernel[i];
        }

        // now normalize each element
        for (int i = 0; i < kernel_size; i++)
        {
            kernel[i] /= sum;
        }

        // convolve with input image horizontally first
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float running_sum = 0.0f;
                for (int i = 0; i < kernel_size; i++)
                {
                    int x_val = x + i - half; // this is the x_val we will take the pixel val from
                    if (x_val < 0)            // fix boundaries
                    {
                        x_val = 0;
                    }
                    if (x_val >= width)
                    {
                        x_val = width - 1;
                    }
                    running_sum += input[y * width + x_val] * kernel[i];
                }
                temp[y * width + x] = running_sum;
            }
        }

        // then vertically
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float running_sum = 0.0f;
                for (int i = 0; i < kernel_size; i++)
                {
                    int y_val = y + i - half; // this is the y_val we will take the pixel val from
                    if (y_val < 0)            // fix boundaries
                    {
                        y_val = 0;
                    }
                    if (y_val >= height)
                    {
                        y_val = height - 1;
                    }
                    running_sum += temp[y_val * width + x] * kernel[i];
                }
                temp2[y * width + x] = running_sum;
            }
        }
        // remap to 0-255 for brightness
        for (int i = 0; i < width * height; i++)
        {
            output[i] = (uint8_t) (255.0f * temp2[i]);
        }
        free(temp);
        free(temp2);
        free(kernel);
        
    }

    // apply airy psf (https://en.wikipedia.org/wiki/Airy_disk)
    else if (psf_type == AIRY)
    {
        float first_min = 1.22f * wavelength / (2.0f * NA); // nm
                  
        // (https://www.keyence.eu/ss/products/microscope/microscope_glossary/terminology/numerical_aperture.jsp)
        // float second_min = 2.233f * wavelength / (2.0f * NA); //nm UNUSED
        // float third_min = 3.238f * wavelength / (2.0f * NA); //nm UNUSED

        int kernel_size = (int) (3 * first_min / pixel_size); // reduce to others if too slow
        if (kernel_size % 2 == 0)
        { 
            kernel_size++; // make sure kernel is odd sized
        }

        int half = kernel_size / 2;

        float *kernel = malloc(kernel_size * kernel_size * sizeof(float));
        float sum = 0.0f;

        float alpha = (2.0f * PI * NA) / wavelength; // nm^-1

        for (int y = 0; y < kernel_size; y++)
        {
            for (int x = 0; x < kernel_size; x++)
            {

                float r = sqrt((x - half) * (x - half) + (y - half) * (y - half)) * pixel_size; // nm

                if (r == 0)
                {
                    kernel[y * kernel_size + x] =
                        1.0f; // fix the singularity that would come from divide by 0
                }
                else
                {
                    float arg = alpha * r; // 2*pi*NA/wavelength * r   
                    kernel[y * kernel_size + x] = pow(
                        (2.0f * j1(arg)) / arg, 2); // airy intensity = (2*j1(x)/x)**2 where j1 is
                                                    // bessel function of first kind (order 1)
                }
                sum += kernel[y * kernel_size + x];
            }
        }
        // normalize
        for (int i = 0; i < kernel_size * kernel_size; i++)
        {
            kernel[i] /= sum;
        }

        // 2D convolution - sadly not separable into x and y
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float running_sum = 0.0f;
                for (int ky = 0; ky < kernel_size; ky++)
                {

                    int y_val = y + ky - half; // this is the y_val we will take the pixel val from
                    if (y_val < 0)
                    {
                        y_val = 0;
                    }
                    if (y_val >= height)
                    {
                        y_val = height - 1;
                    }

                    for (int kx = 0; kx < kernel_size; kx++)
                    {

                        int x_val =
                            x + kx - half; // this is the x_val we will take the pixel val from
                        if (x_val < 0)     // fix boundaries
                        {
                            x_val = 0;
                        }

                        if (x_val >= width)
                        {
                            x_val = width - 1;
                        }

                        running_sum += input[y_val * width + x_val] * kernel[ky * kernel_size + kx];
                    }
                }
                output[y * width + x] = (uint8_t) (255.0f * running_sum);
            }
        }

        free(kernel);
    }
    return;
}
