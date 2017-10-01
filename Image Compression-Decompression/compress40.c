/* compress40.c
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Can compress a ppm/pnm image into a compressed image of codewords. Where
 * the codewords contain the color data for a 2x2 pixel block, and where
 * that codeword is based on the JPEG DCT compression format.
 *
 * Can decompress a file with the codewords next to each other back into
 * a pnm image
 * */
#include <string.h>
#include <stdlib.h>
#include <compress40.h>
#include <assert.h>
#include <except.h>
#include <a2methods.h>
#include <a2blocked.h>
#include <pnm.h>
#include <mem.h>
#include "ppm_manip.h"
#include "pixels.h"
#include "bitpack.h"

#define BLOCKSIZE 2
#define DENOMINATOR 255
#define CHARWIDTH 8
#define NUMCHARIN32BIT 4
#define WORDSIZE 32
/*compression functions*/
static Pnm_ppm trim_image(Pnm_ppm orig_img, A2Methods_T methods);
/*rgb->component, uses cmp_cl*/
static A2Methods_UArray2 conv_to_cmp(Pnm_ppm orig_img,
                                     A2Methods_T methods);
static void rgb_to_cmp_map(int col, int row, A2Methods_UArray2 array,
                    void *elem, void *vcl);
/*component->dct, uses dct_cl*/
static A2Methods_UArray2 cmp_to_dct(A2Methods_UArray2 cmp_array, A2Methods_T methods);
static void cmp_to_dct_map(int col, int row, A2Methods_UArray2 cmp_array,
                    void *elem, void *vcl);
/*dct->codeword*/
static void print_words(A2Methods_UArray2 dct_array, A2Methods_T methods);
static void print_codeword(uint64_t word);

/*decompression functions*/
/*codeword->dct*/
static A2Methods_UArray2 read_in_to_dct(FILE *input, A2Methods_T methods);
/*dct->component, uses inv_dct_cl*/
static A2Methods_UArray2 dct_to_cmp(A2Methods_UArray2 dct_array, A2Methods_T methods);
static void dct_to_cmp_map(int col, int row, A2Methods_UArray2 dct_array,
                    void *elem, void *vcl);
/*component->pnm, uses inv_cmp_cl*/
static Pnm_ppm cmp_to_img(A2Methods_UArray2 cmp_array, A2Methods_T methods);
static void cmp_to_img_map(int col, int row, A2Methods_UArray2 cmp_array,
                    void *elem, void *vcl);

static Pnm_ppm create_ppm(A2Methods_T methods, A2Methods_UArray2 arr);

/*closures for each mapping function*/
struct cmp_cl {
        A2Methods_T methods;
        A2Methods_UArray2 array;
        unsigned denominator;
};

struct dct_cl {
        A2Methods_T methods;
        A2Methods_UArray2 dct_array;
        Cmpnt_Block block;
        int counter;
};

struct inv_dct_cl {
        A2Methods_T methods;
        A2Methods_UArray2 cmp_array;
};

struct inv_cmp_cl {
        A2Methods_T methods;
        A2Methods_UArray2 rgb_array;
};

/* Purpose: compress controller function
 * Does:    1. Reads in image from the file pointer
 *          2. Trims the image so it has an even number of rows and cols
 *          3. Converts the image from Pnm format(rgb) to component
 *             format (Y, Pb, Pr)
 *          4. Converts the component image to a dct image. Where each
 *             element in the dct image represents a 2x2 block in the
 *             component/rgb images
 *          5. Prints the codewords in row major order from the dct image
 */
extern void compress40  (FILE *input)
{
        assert (input != NULL);
        Pnm_ppm orig_img, trim_img;
        A2Methods_UArray2 cmp_img, dct_img;

        orig_img = NULL;
        trim_img = NULL;
        cmp_img  = NULL;
        dct_img  = NULL;

        A2Methods_T methods = uarray2_methods_blocked;

        TRY
                orig_img = Pnm_ppmread(input, methods);
        EXCEPT(Pnm_Badformat)
                fprintf(stderr, "Error: Badly formatted input file.\n");
        END_TRY;

        trim_img = trim_image(orig_img, methods);
        Pnm_ppmfree(&orig_img);

        cmp_img  = conv_to_cmp(trim_img, methods);
        Pnm_ppmfree(&trim_img);

        dct_img  = cmp_to_dct(cmp_img, methods);
        methods->free(&cmp_img);

        print_words(dct_img, methods);
        methods->free(&dct_img);

}

/* Purpose: decompress controller function
 * Does:    1. Reads in the stream of codewords, converts this to dct and
 *             puts the dct blocks into a dct array/image
 *          2. Converts the dct image back to component
 *          3. Converts the component image back to rgb, in Pnm format
 *          4. Prints the decompressed image
 * */
extern void decompress40(FILE *input)
{

        assert(input != NULL);
        A2Methods_T methods = uarray2_methods_blocked;
        A2Methods_UArray2 dct_img, cmp_img;
        Pnm_ppm final_img;

        dct_img = NULL;
        cmp_img = NULL;
        final_img = NULL;

        dct_img = read_in_to_dct(input, methods);

        cmp_img = dct_to_cmp(dct_img, methods);
        methods->free(&dct_img);

        final_img = cmp_to_img(cmp_img, methods);
        methods->free(&cmp_img);

        Pnm_ppmwrite(stdout, final_img);
        Pnm_ppmfree(&final_img);

        /*methods->free(&dct_img);
        methods->free(&cmp_img);
        Pnm_ppmfree(&final_img);*/
}

/* Computes if there is an odd width or height, then trims the image
 * to the necessary width and height using trim_ppm from ppm_manip.h
 * */
static Pnm_ppm trim_image(Pnm_ppm orig_img, A2Methods_T methods)
{
        assert(orig_img != NULL);
        unsigned width, height, trim_amt_w, trim_amt_h;

        width = orig_img->width;
        height = orig_img->height;

        assert (width > 1);
        assert (height > 1);

        trim_amt_w = 0;
        trim_amt_h = 0;
        if ((width % BLOCKSIZE) != 0) {
                trim_amt_w = 1;
        }
        if ((height % BLOCKSIZE) != 0) {
                trim_amt_h = 1;
        }

        return trim_ppm(orig_img, methods, trim_amt_w,  trim_amt_h);
}

/* COMPRESSION FUNCTIONS */

/* Controller function for the conversion from rgb to component
 * packages the original image, methods, and the denominator into a cmp_cl
 * Then runs both arrays through the rgb_to_cmp_map function which does
 * the actual conversion
 * */
static A2Methods_UArray2 conv_to_cmp(Pnm_ppm orig_img,
                                     A2Methods_T methods)
{
        A2Methods_mapfun *map = methods->map_default;
        A2Methods_UArray2 cmp_img;

        cmp_img = methods->new_with_blocksize(orig_img->width,orig_img->height,
                               sizeof(struct Cmpnt_Pixel), BLOCKSIZE);
        struct cmp_cl cl = { methods, cmp_img, orig_img->denominator };
        map(orig_img->pixels, rgb_to_cmp_map, &cl);

        return cmp_img;
}

/* For every pixel in the old rgb image, this converts those pixels to
 * component and stores them in the cmp_img in the closure
 * */
static void rgb_to_cmp_map(int col, int row, A2Methods_UArray2 array,
                    void *elem, void *vcl)
{
        Cmpnt_Pixel *cp;
        struct cmp_cl *cl;
        A2Methods_T m;
        A2Methods_UArray2 cmp_img;
        unsigned denominator;

        (void)array;

        cl      = (struct cmp_cl*) vcl;
        cmp_img = cl->array;
        m       = cl->methods;
        denominator = cl->denominator;

        cp = m->at(cmp_img, col, row);
        *cp = rgb_to_cmpnt((Pnm_rgb)elem, denominator);
}

/* Controller function for the conversion from component to dct blocks
 * packages the cmp image, methods, the counter, and a block  into a dct_cl
 * Then runs both arrays through the cmp_to_dct_map function which does
 * the actual conversion
 * */
static A2Methods_UArray2 cmp_to_dct(A2Methods_UArray2 cmp_array,
                                    A2Methods_T methods)
{
        int new_w, new_h, counter;
        Cmpnt_Block block;

        A2Methods_mapfun *map = methods->map_default;
        A2Methods_UArray2 dct_img;
        new_w = methods->width(cmp_array)  / BLOCKSIZE;
        new_h = methods->height(cmp_array) / BLOCKSIZE;
        dct_img = methods->new(new_w, new_h, sizeof(struct Dct_c));
        counter = 0;

        struct dct_cl cl = { methods, dct_img, block, counter};
        map(cmp_array, cmp_to_dct_map, &cl);
        return dct_img;
}

/* For every pixel in the old cmp image, this converts those pixels to
 * dct and stores them in the dct_img in the closure
 * Note: This stores 4 pixels in a row into the block within the closure
 *       and then converts to a dct on every 4th function call
 * */
static void cmp_to_dct_map(int col, int row, A2Methods_UArray2 cmp_array,
                    void *elem, void *vcl)
{
        double b,c;
        (void)cmp_array;
        struct dct_cl *cl;
        cl = vcl;
        A2Methods_T m;
        m = cl->methods;

        Cmpnt_Pixel *cmp_pix;
        cmp_pix = elem;
        Cmpnt_Pixel pix;
        pix = *cmp_pix;

        Dct_c *dctp;

        cl->counter++;
        if (cl->counter == 1) {
                cl->block.p1 = pix;
        } else if (cl->counter == 2) {
                cl->block.p2 = pix;
        } else if (cl->counter == 3) {
                cl->block.p3 = pix;
        } else if (cl->counter == 4) {
                cl->block.p4 = pix;
                dctp = m->at(cl->dct_array, col / BLOCKSIZE, row / BLOCKSIZE);
                *dctp = dct_y(cl->block);
                /*for some reason, b and c end up switched, we have no idea
                 * why, switching them back here resolves this*/
                b = dctp->b;
                c = dctp->c;
                dctp->b = c;
                dctp->c = b;

                cl->counter = 0;
        }
}

/* Loops through every pixel in the dct_array in row major order, and prints
 * out the codeword using the print_codeword() helper function. Inserts the
 * proper header for the CIF
 * */
static void print_words(A2Methods_UArray2 dct_array, A2Methods_T methods)
{
        /*The mapping function cannot be used because the array can only map
         * by block, while mapping by row is needed
         */
        unsigned width = methods->width(dct_array);
        unsigned height = methods->height(dct_array);
        printf("COMP40 Compressed image format 2\n%u %u", width * BLOCKSIZE,
                                                          height * BLOCKSIZE);
        printf("\n");
        Dct_c *dct_pix;
        Dct_c d;
        uint64_t word;

        for (unsigned j = 0; j < height; j++) {
                /*col varies more rapidly*/
                for (unsigned i = 0; i < width; i++) {
                        dct_pix = methods->at(dct_array, i, j);
                        d = *dct_pix;
                        word = create_codeword(d);
                        print_codeword(word);
                }
        }
}

/*converts every 8 bits into a char to be printed to the file. Does so in
 * big endian order
 * */
static void print_codeword(uint64_t word)
{
        for (int i = WORDSIZE - CHARWIDTH; i >= 0; i -= CHARWIDTH) {
                putchar((int) Bitpack_getu(word, CHARWIDTH, i));
        }
}

/*DECOMPRESSION FUNCTOINS */

/* Reasd in the header file, makes a new dct_array, and puts in the extracted
 * dct from every codeword into the new dct_array in row major order
 * Notes: reads in a character in big endian format. 4 characters are read in
 *        make up a 32 bit codeword
 * */
static A2Methods_UArray2 read_in_to_dct(FILE *input, A2Methods_T methods)
{
        unsigned width, height;
        int read = fscanf(input, "COMP40 Compressed image format 2\n%u %u"
                               , &width, &height);
        width = width / BLOCKSIZE;
        height = height / BLOCKSIZE;
        assert(read == 2);
        int c = getc(input);
        assert(c == '\n');
        int cur_char, lsb;
        uint64_t cur_word = 0;
        Dct_c *dct_pix;
        Dct_c pix;

        A2Methods_UArray2 dct_array = methods->new((int)width,(int) height,
                                                   sizeof(struct Dct_c));
        for (unsigned j = 0; j < height; j++) {
                /*col varies more rapidly*/
                for (unsigned i = 0; i < width; i++) {
                        cur_word = 0;
                        for (int k = 0; k < NUMCHARIN32BIT; k++) {
                                cur_char = fgetc(input);
                                lsb = WORDSIZE - ((k + 1) * CHARWIDTH);
                                cur_word = Bitpack_newu(cur_word, CHARWIDTH,
                                                        lsb, cur_char);
                        }
                        pix = extract_codeword(cur_word);
                        dct_pix = methods->at(dct_array, i, j);
                        *dct_pix = pix;
                }
        }
        return dct_array;
}

/* Controller function for the conversion from dct blocks to component
 * packages the cmp image and the A2Methods into a inv_dct_cl
 * Then runs both arrays through the dct_to_cmp_map function which does
 * the actual conversion
 * */
static A2Methods_UArray2 dct_to_cmp(A2Methods_UArray2 dct_array,
                                    A2Methods_T methods)
{
        int new_w, new_h;
        A2Methods_mapfun *map = methods->map_default;
        A2Methods_UArray2 cmp_array;
        new_w = methods->width(dct_array) * BLOCKSIZE;
        new_h = methods->height(dct_array) * BLOCKSIZE;
        cmp_array = methods->new(new_w, new_h, sizeof(struct Cmpnt_Pixel));

        struct inv_dct_cl cl = { methods, cmp_array};
        map(dct_array, dct_to_cmp_map, &cl);
        return cmp_array;
}

/* For every pixel in the old dct array, this converts those blocks to
 * 4 component pixels and stores them in the component array in the closure
 * Note: This stores 4 pixels in a block in the component array every time
 * */
static void dct_to_cmp_map(int col, int row, A2Methods_UArray2 dct_array,
                    void *elem, void *vcl)
{
        (void)dct_array;
        struct inv_dct_cl *cl;
        A2Methods_T m;
        A2Methods_UArray2 cmp_img;
        Dct_c *dct_pix;
        Dct_c pix;
        Cmpnt_Block block;
        Cmpnt_Pixel *cp;

        cl = vcl;
        m  = cl->methods;
        cmp_img = cl->cmp_array;
        dct_pix = elem;
        pix = *dct_pix;
        block = inverse_dct(pix);

        int start_col = col * BLOCKSIZE;
        int start_row = row * BLOCKSIZE;

        cp = m->at(cmp_img, start_col, start_row);
        *cp = block.p1;
        cp = m->at(cmp_img, start_col + 1, start_row);
        *cp = block.p2;
        cp = m->at(cmp_img, start_col, start_row + 1);
        *cp = block.p3;
        cp = m->at(cmp_img, start_col + 1, start_row + 1);
        *cp = block.p4;

}

/* Controller function for the conversion from component to an rgb Ppm image
 * packages the image, and A2methods into a inv_comp_cl
 * Then runs both arrays through the cmp_to_img_map function which does
 * the actual conversion
 * Returns a Pnm_ppm image using the helper function create_ppm)
 * */
static Pnm_ppm cmp_to_img(A2Methods_UArray2 cmp_array, A2Methods_T methods)
{
        int w = methods->width(cmp_array);
        int h = methods->height(cmp_array);
        A2Methods_mapfun *map = methods->map_default;
        A2Methods_UArray2 rgb_array;

        rgb_array = methods->new(w, h, sizeof(struct Pnm_rgb));

        struct inv_cmp_cl cl = { methods, rgb_array };
        map(cmp_array, cmp_to_img_map, &cl);
        return create_ppm(methods, rgb_array);
}

/* For every pixel in the old cmp array, this converts those  to
 * rgb pixels and stores them in the rgb array in the closure
 * */
static void cmp_to_img_map(int col, int row, A2Methods_UArray2 cmp_array,
                    void *elem, void *vcl)
{
        Pnm_rgb p;
        Cmpnt_Pixel *cp;
        struct inv_cmp_cl *cl;
        A2Methods_T m;
        A2Methods_UArray2 rgb_array;
        (void)cmp_array;

        cl      = vcl;
        rgb_array = cl->rgb_array;
        m       = cl->methods;

        cp = elem;
        p = m->at(rgb_array, col, row);
        *p = cmpnt_to_rgb(*cp, DENOMINATOR);
}

static Pnm_ppm create_ppm(A2Methods_T methods, A2Methods_UArray2 arr)
{
        Pnm_ppm new_ppm;
        NEW(new_ppm);
        new_ppm->width = methods->width(arr);
        new_ppm->height = methods->height(arr);
        new_ppm->pixels = arr;
        new_ppm->denominator = DENOMINATOR;
        new_ppm->methods = methods;

        return new_ppm;
}

#undef BLOCKSIZE
#undef DENOMINATOR
#undef CHARWIDTH
#undef NUMCHARIN32BIT
#undef WORDSIZE

