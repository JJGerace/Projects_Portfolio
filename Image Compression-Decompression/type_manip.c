/* type_manip.c
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * type_manip.c converts the 4 parameters a,b,c,d to and from
 * their relevant integer type, can also be reused to work with
 * other conversions (ex: a double to a 10 bit integer instead of
 * a 9 bit) by changing constants
 * */
#include <stdint.h>
#define ASCALE 63
#define BCDSCALE 100
/*.5 is added to make casting as an integer round to the nearest integer,
 * this prevents always having floor rounding
 * */
int64_t double_to_signed(double d)
{
        return (int64_t)(d + 0.5);
}

/*abcd are rounded to 0.3 or -0.3 if too large, then converted to a range of
 * -15 to 15*/
int64_t double_bcd_to_signed(double d)
{
        if (d < -0.3) {
                d = -0.3;
        } else if (d > 0.3) {
                d = 0.3;
        }

        int num = (int)((d * BCDSCALE) + .5);
        return num;
}

/*a is forced to be non negative, then converted to a scaled unsinged
 * integer
 * */
uint64_t double_a_to_unsigned(double d)
{
        if (d < 0)
                d = 0;
        uint64_t num = (uint64_t)((d * ASCALE) + 0.5);
        return num;
}

double signed_bcd_to_double(int64_t i)
{
        double d;
        d = (double)i;
        d = (d - .5) / BCDSCALE;
        return d;
}

double unsigned_a_to_double(uint64_t i)
{
        double d;
        d = (double)i;
        d = (d - .5) / ASCALE;
        return d;
}
#undef ASCALE
#undef BCDSCALE

