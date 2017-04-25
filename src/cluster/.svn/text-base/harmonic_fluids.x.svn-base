/*
 * harmonic_fluids.x
 *
 * SUN RPC interface for harmonic fluids computation
 *
 * Changxi Zheng (2008)
 */

struct position
{
    double x;
    double y;
    double z;
};

typedef struct position* position_ptr;

/****************** boundary samples ****************/
struct hf_samples
{
    int             free_len;
    int             left_len;
    int             right_len;
    int             front_len;
    int             back_len;
    int             bottom_len;

    position_ptr    free_sample;
    position_ptr    free_normal;
    position_ptr    left_sample;
    position_ptr    right_sample;
    position_ptr    front_sample;
    position_ptr    back_sample;
    position_ptr    bottom_sample;
};

/****************** bubble info for multipole estimation ***************/
struct hf_bubble
{
    int             id;     /* bubble id */
    double          freq;   /* frequency */
    double          rad;    /* bubble radius */
    position_ptr    pos;
    position_ptr    regular_src_in<>;
    position_ptr    regular_src_out<>;
};

typedef hf_bubble hf_bubbles<>;

/******************** Results *******************/
struct complex
{
    double real;
    double imag;
};

struct hf_result
{
    int     id;     /* bubble id */
    complex c<>;
};

union hf_results switch (int errno) 
{
    case 0:
        hf_result data<>;
    default:
        void;
};
/******************** Service interface ******************/
program HARMONIC_FLUIDS_PROGRAM
{
    version HARMONIC_FLUIDS_VERSION_1
    {
        int     TEST(hf_bubbles) = 0;
        /*
        int     HF_CONFIGURE(hf_samples) = 1;
        */
        hf_results HF_SOLVE(hf_bubbles) = 2;
    } = 1;
} = 0x20000586;
