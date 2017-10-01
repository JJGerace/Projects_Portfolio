/* pixels.h
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Provides conversions between rgb pixels, component pixels, dct pixel blocks
 * and the codeword using bitpacker
 * */
#include <pnm.h>
#include <bitpack.h>
/*structs represent the intermediate pixels */
/* component pixels */
typedef struct Cmpnt_Pixel{
        double Y;
        double Pb;
        double Pr;
} Cmpnt_Pixel;

/* discrete cosine transformed block of pixels */
typedef struct Dct_c {
        double a;
        double b;
        double c;
        double d;
        double _Pb;
        double _Pr;
} Dct_c;

typedef struct Cmpnt_Block {
        Cmpnt_Pixel p1, p2, p3, p4;
} Cmpnt_Block;
/*conversion functions*/
Cmpnt_Pixel rgb_to_cmpnt(Pnm_rgb rgb, unsigned denominator);
struct Pnm_rgb cmpnt_to_rgb(Cmpnt_Pixel cmpnt, unsigned denominator);

Dct_c dct_y(Cmpnt_Block pblock);
Cmpnt_Block inverse_dct(Dct_c dct);

/*helper conversion functions */
double get_Pb_avg(double pb1, double pb2, double pb3, double pb4);
double get_Pr_avg(double pr1, double pr2, double pr3, double pr4);

/*conversion to codeword functions*/
uint64_t create_codeword(Dct_c components);
Dct_c extract_codeword(uint64_t word);
