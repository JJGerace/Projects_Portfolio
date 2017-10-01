/* ppm_manip.c
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Set up to provide a variety of ppm manipualtions. Though only a trim
 * function has been implemented. Which trims a givem ppm image by a
 * given amount in each dimension
 * */
#include <pnm.h>
#include <mem.h>
#include <a2plain.h>
#include <assert.h>
#include "ppm_manip.h"



Pnm_ppm trim_ppm(Pnm_ppm img, A2Methods_T m,
                  unsigned trim_amt_w, unsigned trim_amt_h)
{
        assert (img != NULL);
        assert (m != NULL);
        assert (img->width > trim_amt_w);
        assert (img->height > trim_amt_h);
        unsigned new_w = img->width - trim_amt_w;
        unsigned new_h = img->height - trim_amt_h;

        Pnm_ppm newppm;
        NEW(newppm);
        unsigned i, j;
        A2Methods_UArray2 newarr;
        newarr = m->new(new_w, new_h, sizeof(struct Pnm_rgb));
        /*ignores edge pixels the client wants to trim if trim_amt_h and/or
         * trim_amt_w are nonzero
         * */
        for (i = 0; i < new_w ; ++i) {
                for (j = 0; j < new_h; ++j) {
                        Pnm_rgb p, q;
                        p = m->at(img->pixels, i, j);
                        q = m->at(newarr, i, j);
                        *q = *(Pnm_rgb)p;
                }
        }
        newppm->width = new_w;
        newppm->height = new_h;
        newppm->denominator = img->denominator;
        newppm->pixels = newarr;
        newppm->methods = m;

        return newppm;
}
