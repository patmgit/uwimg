#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

#define CLAMP(upper, lower, n) ((n > upper) ? upper : ((n < lower) ? lower : n))

// helper method to calculate the pixel index in the data array
int calculate_pixel(int x, int w, int y, int h, int c)
{
    return x + (y*w) + (c*w*h);
}

// helper method to calculate the grayscale value
float luma_helper(image im, int x, int y)
{
    return (0.299 * get_pixel(im, x, y, 0)) + (0.587 * get_pixel(im, x, y, 1)) + (0.114 * get_pixel(im, x, y, 2));
}

// helper method to change every pixel in a channel with the passed in function
void set_helper(image src, image dest, int c, float (*f)(image, int, int))
{
    int x, y;
    for (y = 0; y < dest.h; y++) {
      for (x = 0; x < dest.w; x++) {
        set_pixel(dest, x, y, c, f(src, x, y));
      }
    }
}

float get_pixel(image im, int x, int y, int c)
{
    int clampedx = CLAMP(im.w - 1, 0, x);
    int clampedy = CLAMP(im.h - 1, 0, y);
    int clampedc = CLAMP(im.c - 1, 0, c);
    return im.data[calculate_pixel(clampedx, im.w, clampedy, im.h, clampedc)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (c > im.c || y > im.h || x > im.w) {
      return;
    } else {
      im.data[calculate_pixel(x, im.w, y, im.h, c)] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, sizeof(float) * (im.w * im.h * im.c));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    set_helper(im, gray, 0, luma_helper);
    return gray;
}

void shift_image(image im, int c, float v)
{
    // nested helper function that returns a shifted pixel value
    float shift_helper(image img, int x, int y) { return get_pixel(img, x, y, c) + v; }
    set_helper(im, im, c, shift_helper);
}

void clamp_image(image im)
{
    // nested nested helper function that clamps the pixel
    void clamp_helper(int c) {
      float clamp_helper_helper(image img, int x, int y)
      { 
        return CLAMP(1, 0, get_pixel(img, x, y, c));
      }
      set_helper(im, im, c, clamp_helper_helper);
    }
    clamp_helper(0);
    clamp_helper(1);
    clamp_helper(2);
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{ 
    int x, y;
    float r, g, b, v, m, c, s, h, hprime;
    for (y = 0; y < im.h; y++) {
      for (x = 0; x < im.w; x++) {
        r = get_pixel(im, x, y, 0);
        g = get_pixel(im, x, y, 1);
        b = get_pixel(im, x, y, 2);
        v = three_way_max(r, g, b);
        m = three_way_min(r, g, b);
        c = v - m;
        s = ((v == 0) ? 0 : (c / v));
        if (c == 0) {
          h = 0;
        } else {
          if (v == r) {
            hprime = (g - b) / c;
          } else if (v == g) {
            hprime = ((b - r) / c) + 2;
          } else {
            hprime = ((r - g) / c) + 4;
          }
          h = hprime / 6;
          if (hprime < 0) {
            h = h + 1;
          }
        }
        set_pixel(im, x, y, 0, h);
        set_pixel(im, x, y, 1, s); 
        set_pixel(im, x, y, 2, v);
      }
    }
}

void hsv_to_rgb(image im)
{
    int x, y;
    float h, s, v, c, i, m, hprime, r, g, b;
    for (y = 0; y < im.h; y++) {
      for (x = 0; x < im.w; x++) {
        h = get_pixel(im, x, y, 0);
        s = get_pixel(im, x, y, 1);
        v = get_pixel(im, x, y, 2); 
        c = v * s;
        hprime = h * 6.0f;
        i = c * (1.0f - fabs(fmod(hprime, 2.0f) - 1.0f));
        m = v - c;
        if (hprime <= 1) {
          r = c;
          g = i;
          b = 0;
        } else if (hprime <= 2) {
          r = i;
          g = c;
          b = 0;
        } else if (hprime <= 3) {
          r = 0;
          g = c;
          b = i;
        } else if (hprime <= 4) {
          r = 0;
          g = i;
          b = c;
        } else if (hprime <= 5) {
          r = i;
          g = 0;
          b = c;
        } else {
          r = c;
          g = 0;
          b = i;
        }
        r = r + m;
        g = g + m;
        b = b + m;       
        set_pixel(im, x, y, 0, r);
        set_pixel(im, x, y, 1, g); 
        set_pixel(im, x, y, 2, b);
      }
    }
}

void scale_image(image im, int c, float v)
{
    // nested helper function that returns a scaled pixel value
    float scale_helper(image img, int x, int y) { return get_pixel(img, x, y, c) * v; }
    set_helper(im, im, c, scale_helper);
}
