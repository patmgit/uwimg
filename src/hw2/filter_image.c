#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float sum = 0.0f;
    for (int c = 0; c < im.c; c++) {
      for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
          sum += get_pixel(im, x, y, c);
        }
      }
    }

    for (int c = 0; c < im.c; c++) {
      for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
          set_pixel(im, x, y, c, (get_pixel(im, x, y, c) / sum));
        }
      }
    }
}

image make_box_filter(int w)
{
    image im = make_image(w, w, 1);

    for (int y = 0; y < im.h; y++) {
      for (int x = 0; x < im.w; x++) {
        set_pixel(im, x, y, 0, 1.0f);
      }
    }
    l1_normalize(im);

    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // assume filter is odd?
    
    assert(filter.c == im.c || filter.c == 1);
    image res;
    if (preserve == 1) {
       res = make_image(im.w, im.h, im.c);
    } else {
       res = make_image(im.w, im.h, 1);
    }
 
    int center_x = filter.w / 2;
    int center_y = filter.h / 2;

    for (int x = 0; x < res.w; x++) {
      for (int y = 0; y < res.h; y++) {
        for (int c = 0; c < res.c; c++) {
          float value = 0;
          for (int filter_x = x - center_x; filter_x <= x + center_x; filter_x++) {
            for (int filter_y = y - center_y; filter_y <= y + center_y; filter_y++) {
              if (preserve == 1) {
                if (filter.c == 1) {
                  // filter has 1, im has 3, res has 3
                  value += get_pixel(im, filter_x, filter_y, c) * get_pixel(filter, (filter_x - (x - center_x)), (filter_y - (y - center_y)), 0);
                } else {
                  // filter has 3, im has 3, res has 3
                  value += get_pixel(im, filter_x, filter_y, c) * get_pixel(filter, (filter_x - (x - center_x)), (filter_y - (y - center_y)), c);
                }
              } else {
                for (int image_c = 0; image_c < im.c; image_c++) {
                  if (filter.c == 1) {
                    // filter has 1, im has 3, res has 1
                    value += get_pixel(im, filter_x, filter_y, image_c) * get_pixel(filter, (filter_x - (x - center_x)), (filter_y - (y - center_y)), c);
                  } else {
                    // filter has 3, im has 3, res has 1
                    value += get_pixel(im, filter_x, filter_y, image_c) * get_pixel(filter, (filter_x - (x - center_x)), (filter_y - (y - center_y)), image_c);
                  }
                }
              }
            }
          }
          set_pixel(res, x, y, c, value);
        }
      }
    }
    
    return res;
}

image make_highpass_filter()
{
    image im = make_image(3, 3, 1);
    
    set_pixel(im, 0, 0, 0, 0.0f);
    set_pixel(im, 1, 0, 0, -1.0f);
    set_pixel(im, 2, 0, 0, 0.0f);
    set_pixel(im, 0, 1, 0, -1.0f);
    set_pixel(im, 1, 1, 0, 4.0f);
    set_pixel(im, 2, 1, 0, -1.0f);
    set_pixel(im, 0, 2, 0, 0.0f);
    set_pixel(im, 1, 2, 0, -1.0f);
    set_pixel(im, 2, 2, 0, 0.0f);

    return im;
}

image make_sharpen_filter()
{
    image im = make_image(3, 3, 1);
    
    set_pixel(im, 0, 0, 0, 0.0f);
    set_pixel(im, 1, 0, 0, -1.0f);
    set_pixel(im, 2, 0, 0, 0.0f);
    set_pixel(im, 0, 1, 0, -1.0f);
    set_pixel(im, 1, 1, 0, 5.0f);
    set_pixel(im, 2, 1, 0, -1.0f);
    set_pixel(im, 0, 2, 0, 0.0f);
    set_pixel(im, 1, 2, 0, -1.0f);
    set_pixel(im, 2, 2, 0, 0.0f);

    return im;
}

image make_emboss_filter()
{
    image im = make_image(3, 3, 1);
    
    set_pixel(im, 0, 0, 0, -2.0f);
    set_pixel(im, 1, 0, 0, -1.0f);
    set_pixel(im, 2, 0, 0, 0.0f);
    set_pixel(im, 0, 1, 0, -1.0f);
    set_pixel(im, 1, 1, 0, 1.0f);
    set_pixel(im, 2, 1, 0, 1.0f);
    set_pixel(im, 0, 2, 0, 0.0f);
    set_pixel(im, 1, 2, 0, 1.0f);
    set_pixel(im, 2, 2, 0, 2.0f);

    return im;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should use preserve on sharpen and emboss and not for highpass because the purpose of highpass is to find edges, so we just want a one channel image. On the other hand, sharpen and emboss are effects we want to apply to the image itself, so we should preserve all three channels.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Because I was not able to make and test, my guess is that sharpen and emboss need to be clamped because both filters have matrices which sum to 1, meaning that the output may exceed the highest value.

image make_gaussian_filter(float sigma)
{
    int size = (int) ceilf(sigma * 6);
    if (size % 2 == 0) {
      size = size + 1;
    }
    image im = make_image(size, size, 1);

    for (int x = - im.w / 2; x <= im.w / 2; x++) {
      for (int y = - im.w / 2; y <= im.w / 2; y++) {
        float term = exp((-1 * (pow(x, 2) + pow (y, 2))) / (2 * pow(sigma, 2)));
        term /= (TWOPI * pow(sigma, 2));
        set_pixel(im, x + (im.w / 2), y + (im.h / 2), 0, term);
      }
    }
    l1_normalize(im);
    
    return im;
}

image add_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    
    image im = make_image(a.w, a.h, a.c);
   
    for (int x = 0; x < a.w; x++) {
      for (int y = 0; y < a.h; y++) {
        for (int c = 0; c < a.c; c++) {
          set_pixel(im, x, y, c, (get_pixel(a, x, y, c) + get_pixel(b, x, y, c)));
        }
      }
    }
    
    return im;
}

image sub_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    
    image im = make_image(a.w, a.h, a.c);
    
    for (int x = 0; x < a.w; x++) {
      for (int y = 0; y < a.h; y++) {
        for (int c = 0; c < a.c; c++) {
          set_pixel(im, x, y, c, (get_pixel(a, x, y, c) - get_pixel(b, x, y, c)));
        }
      }
    }
    
    return im;
}

image make_gx_filter()
{
    image im = make_image(3, 3, 1);
    
    set_pixel(im, 0, 0, 0, -1.0f);
    set_pixel(im, 1, 0, 0, 0.0f);
    set_pixel(im, 2, 0, 0, 1.0f);
    set_pixel(im, 0, 1, 0, -2.0f);
    set_pixel(im, 1, 1, 0, 0.0f);
    set_pixel(im, 2, 1, 0, 2.0f);
    set_pixel(im, 0, 2, 0, -1.0f);
    set_pixel(im, 1, 2, 0, 0.0f);
    set_pixel(im, 2, 2, 0, 1.0f);

    return im;
}

image make_gy_filter()
{
    image im = make_image(3, 3, 1);
    
    set_pixel(im, 0, 0, 0, -1.0f);
    set_pixel(im, 1, 0, 0, -2.0f);
    set_pixel(im, 2, 0, 0, -1.0f);
    set_pixel(im, 0, 1, 0, 0.0f);
    set_pixel(im, 1, 1, 0, 0.0f);
    set_pixel(im, 2, 1, 0, 0.0f);
    set_pixel(im, 0, 2, 0, 1.0f);
    set_pixel(im, 1, 2, 0, 2.0f);
    set_pixel(im, 2, 2, 0, 1.0f);

    return im;
}

void feature_normalize(image im)
{
    float min = get_pixel(im, 0, 0, 0);
    float max = get_pixel(im, 0, 0, 0);
    float current = get_pixel(im, 0, 0, 0);
    for (int x = 0; x < im.w; x++) {
      for (int y = 0; y < im.h; y++) {
        for (int c = 0; c < im.c; c++) {
          current = get_pixel(im, x, y, c);
          if (current < min) {
            min = current;
          }
          if (current > max) {
            max = current;
          }
        }
      }
    }
    float range = max - min;
    if (range == 0) {
      for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
          for (int c = 0; c < im.c; c++) {
            set_pixel(im, x, y, c, 0.0f);
          }
        }
      }
    } else {
      for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
          for (int c = 0; c < im.c; c++) {
            current = ((get_pixel(im, x, y, c) - min) / range);
            set_pixel(im, x, y, c, current);
          }
        }
      }
    }
}

image *sobel_image(image im)
{
    image *res = calloc(2, sizeof(image));

    res[0] = make_image(im.w, im.h, 1);
    res[1] = make_image(im.w, im.h, 1);
    
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image gx_im = convolve_image(im, gx_filter, 0);    
    image gy_im = convolve_image(im, gy_filter, 0);

    for (int x = 0; x < im.w; x++) {
      for (int y = 0; y < im.h; y++) {
	float gx = get_pixel(gx_im, x, y, 0);
        float gy = get_pixel(gy_im, x, y, 0);
        float mag = sqrtf(pow(gx, 2) + pow(gy, 2));
        float dir = atan2f(gy, gx);
        set_pixel(res[0], x, y, 0, mag);
        set_pixel(res[1], x, y, 0, dir);
      }
    }

    free_image(gx_filter);
    free_image(gy_filter);
    free_image(gx_im);
    free_image(gy_im);
 
    return res;
}

image colorize_sobel(image im)
{

    image *res = sobel_image(im);
    image mag = res[0];
    image dir = res[1];
    
    image color = make_image(im.w, im.h, 3);

    for (int x = 0; x < im.w; x++) {
      for (int y = 0; y < im.h; y++) {
        set_pixel(color, x, y, 0, get_pixel(mag, x, y, 0));
        set_pixel(color, x, y, 1, get_pixel(mag, x, y, 0));
        set_pixel(color, x, y, 2, get_pixel(dir, x, y, 0));
      }
    }
    hsv_to_rgb(color);

    free_image(mag);
    free_image(dir);
    free(res);

    return color;
}
