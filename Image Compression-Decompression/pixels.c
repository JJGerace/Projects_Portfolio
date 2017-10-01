/* pixels.c
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Provides conversions between rgb pixels, component pixels, dct pixel blocks
 * and the codeword using bitpacker
 * */
#include <pixels.h>
#include <mem.h>
#include <arith40.h>
#include "bitpack.h"
#include "type_manip.h"
#include "a2blocked.h"

#define AWIDTH 6
#define BCDWIDTH 6
#define PBPRWIDTH 4
#define ALSB 26
#define BLSB 20
#define CLSB 14
#define DLSB 8
#define PBLSB 4
#define PRLSB 0

Cmpnt_Pixel rgb_to_cmpnt(Pnm_rgb rgb, unsigned denominator)
{
        double y, pb, pr, r, g, b;
        struct Cmpnt_Pixel cpix;

        r = (double) rgb->red / (double) denominator;
        g = (double) rgb->green / (double) denominator;
        b = (double) rgb->blue / (double) denominator;

        y = 0.299 * r + 0.587 * g + 0.114 * b;
        pb = -0.168736 * r - 0.331264 * g + 0.5 * b;
        pr = 0.5 * r - 0.418688 * g - 0.081312 * b;

        cpix.Y = y;
        cpix.Pb = pb;
        cpix.Pr = pr;

        return cpix;
}

/* The converted rgb values are invalid if less than 0 or greater than the
 * denominator, and are rounded to 0 or denominator in that case
 * */
struct Pnm_rgb cmpnt_to_rgb(Cmpnt_Pixel cmpnt, unsigned denominator)
{
        double r, g, b;
        struct Pnm_rgb rgb;
        r = 1.0 * cmpnt.Y + 0.0 * cmpnt.Pb + 1.402 * cmpnt.Pr;
        g = 1.0 * cmpnt.Y - 0.344136 * cmpnt.Pb - 0.714136 * cmpnt.Pr;
        b = 1.0 * cmpnt.Y + 1.772 * cmpnt.Pb + 0.0 * cmpnt.Pr;
        r = r * denominator;
        g = g * denominator;
        b = b * denominator;

        if (r < 0) {
                r = 0;
        } else if (r > denominator) {
                r = denominator;
        }
        if (g < 0) {
                g = 0;
        } else if (g > denominator) {
                g = denominator;
        }
        if (b < 0) {
                b = 0;
        } else if (b > denominator) {
                b = denominator;
        }

        rgb.red = double_to_signed(r);
        rgb.green = double_to_signed(g);
        rgb.blue = double_to_signed(b);
        return rgb;
}

/* converts a block of component pixels into a dct representation of the block
 * */
Dct_c dct_y(Cmpnt_Block pblock)
{
        double Y1, Y2, Y3, Y4, _pb, _pr, a, b, c, d;
        Dct_c dct;

        Y1 = pblock.p1.Y;
        Y2 = pblock.p2.Y;
        Y3 = pblock.p3.Y;
        Y4 = pblock.p4.Y;

        a = (Y4 + Y3 + Y2 + Y1) / 4.000;
        b = (Y4 + Y3 - Y2 - Y1) / 4.000;
        c = (Y4 - Y3 + Y2 - Y1) / 4.000;
        d = (Y4 - Y3 - Y2 + Y1) / 4.000;

        _pb = get_Pb_avg(pblock.p1.Pb, pblock.p2.Pb,
                         pblock.p3.Pb, pblock.p4.Pb);
        _pr = get_Pr_avg(pblock.p1.Pr, pblock.p2.Pr,
                         pblock.p3.Pr, pblock.p4.Pr);

        dct.a = a;
        dct.b = b;
        dct.c = c;
        dct.d = d;
        dct._Pb = _pb;
        dct._Pr = _pr;

        return dct;
}

/* takes in a dct and returns a component block with inverse dct values */
Cmpnt_Block inverse_dct(Dct_c dct)
{
        double Y1, Y2, Y3, Y4;
        Cmpnt_Block cblock;
        Cmpnt_Pixel p1, p2, p3, p4;

        Y1 = dct.a - dct.b - dct.c + dct.d;
        Y2 = dct.a - dct.b + dct.c - dct.d;
        Y3 = dct.a + dct.b - dct.c - dct.d;
        Y4 = dct.a + dct.b + dct.c + dct.d;

        p1.Y = Y1;
        p2.Y = Y2;
        p3.Y = Y3;
        p4.Y = Y4;

        p1.Pb = dct._Pb;
        p2.Pb = dct._Pb;
        p3.Pb = dct._Pb;
        p4.Pb = dct._Pb;

        p1.Pr = dct._Pr;
        p2.Pr = dct._Pr;
        p3.Pr = dct._Pr;
        p4.Pr = dct._Pr;

        cblock.p1 = p1;
        cblock.p2 = p2;
        cblock.p3 = p3;
        cblock.p4 = p4;

        return cblock;
}

double get_Pb_avg(double pb1, double pb2, double pb3, double pb4)
{
        return ((double) (pb1 + pb2 + pb3 + pb4) / 4.000);
}

double get_Pr_avg(double pr1, double pr2, double pr3, double pr4)
{
        return ((double) (pr1 + pr2 + pr3 + pr4) / 4.000);
}

uint64_t create_codeword(Dct_c components)
{
        uint64_t a_to_pack;
        int64_t b_to_pack, c_to_pack, d_to_pack;
        uint64_t quant_Pb, quant_Pr;
        /*helper functions from type_manip.c and arith40*/
        a_to_pack = double_a_to_unsigned(components.a);
        b_to_pack = double_bcd_to_signed(components.b);
        c_to_pack = double_bcd_to_signed(components.c);
        d_to_pack = double_bcd_to_signed(components.d);
        quant_Pb = Arith40_index_of_chroma((float)components._Pb);
        quant_Pr = Arith40_index_of_chroma((float)components._Pr);

        /*create the word with bitpacker.c*/
        uint64_t word = 0;
        word = Bitpack_newu(word, AWIDTH,    ALSB, a_to_pack);
        word = Bitpack_news(word, BCDWIDTH,  BLSB, b_to_pack);
        word = Bitpack_news(word, BCDWIDTH,  CLSB, c_to_pack);
        word = Bitpack_news(word, BCDWIDTH,  DLSB, d_to_pack);
        word = Bitpack_newu(word, PBPRWIDTH, PBLSB, quant_Pb);
        word = Bitpack_newu(word, PBPRWIDTH, PRLSB, quant_Pr);
        return word;
}

Dct_c extract_codeword(uint64_t word)
{
        uint64_t unpacked_a;
        int64_t unpacked_b, unpacked_c, unpacked_d;
        uint64_t unpacked_Pb, unpacked_Pr;
        /* from bitpacker.c*/
        unpacked_a  = Bitpack_getu(word, AWIDTH, ALSB);
        unpacked_b  = Bitpack_gets(word, BCDWIDTH, BLSB);
        unpacked_c  = Bitpack_gets(word, BCDWIDTH, CLSB);
        unpacked_d  = Bitpack_gets(word, BCDWIDTH, DLSB);
        unpacked_Pb = Bitpack_getu(word, PBPRWIDTH, PBLSB);
        unpacked_Pr = Bitpack_getu(word, PBPRWIDTH, PRLSB);
        /* helper functions from type_manip.c and arith40*/
        Dct_c new_dct;
        new_dct.a   = unsigned_a_to_double(unpacked_a);
        new_dct.b   = signed_bcd_to_double(unpacked_b);
        new_dct.c   = signed_bcd_to_double(unpacked_c);
        new_dct.d   = signed_bcd_to_double(unpacked_d);
        new_dct._Pb = (double)Arith40_chroma_of_index(unpacked_Pb);
        new_dct._Pr = (double)Arith40_chroma_of_index(unpacked_Pr);

        return new_dct;
}

#undef DENOMINATOR
#undef AWIDTH
#undef BCDWIDTH
#undef PBPRWIDTH
#undef ALSB
#undef BLSB
#undef CLSB
#undef DLSB
#undef PBLSB
#undef PRLSB


