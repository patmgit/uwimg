#include <math.h>
#include "image.h"

image resize_helper(image im, int w, int h, float (*f)(image, float, float, int))
{
    float a_x = ((float) im.w)/((float) w);
    float a_y = ((float) im.h)/((float) h);
    float b_x = -0.5f - (-0.5f * a_x);
    float b_y = -0.5f - (-0.5f * a_y); 
    image res = make_image(w, h, im.c);
    int x, y, ch;
    for (ch = 0; ch < res.c; ch++) {
      for (y = 0; y < res.h; y++) {
        for (x = 0; x < res.w; x++) {
          float mapped_x = (a_x * x) + b_x;
          float mapped_y = (a_y * y) + b_y;
          set_pixel(res, x, y, ch, f(im, mapped_x, mapped_y, ch));
        }
      }
    }
    return res;
} 

float nn_interpolate(image im, float x, float y, int c)
{
    return get_pixel(im, (int)round(x), (int)round(y), c);
}

image nn_resize(image im, int w, int h)
{
    return resize_helper(im, w, h, nn_interpolate);
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    float v1 = get_pixel(im, (int)floor(x), (int)floor(y), c);
    float v2 = get_pixel(im, (int)ceil(x), (int)floor(y), c);
    float v3 = get_pixel(im, (int)floor(x), (int)ceil(y), c);
    float v4 = get_pixel(im, (int)ceil(x), (int)ceil(y), c);

    float d1 = x - floor(x);
    float d2 = ceil(x) - x;
    float d3 = y - floor(y);
    float d4 = ceil(y) - y;

    float q1 = (v1 * d2) + (v2 * d1);
    float q2 = (v3 * d2) + (v4 * d1);
    return (q1 * d4) + (q2 * d3);
}

image bilinear_resize(image im, int w, int h)
{
    return resize_helper(im, w, h, bilinear_interpolate);
}

