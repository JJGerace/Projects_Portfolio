/* ppm_manip.h
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Set up to provide a variety of ppm manipualtions. Though only a trim
 * function has been implemented. Which trims a givem ppm image by a
 * given amount in each dimension
 * */
#include <pnm.h>
#include <a2plain.h>

/* trims a given ppm by a given amount in each direction */
Pnm_ppm trim_ppm(Pnm_ppm img, A2Methods_T m,
                  unsigned trim_amt_w, unsigned trim_amt_h);

