
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define X 423
#define Y 523
#define WIDTH 1280
#define HEIGHT 720
#define B 2
#define KERN_SIZE 3
#define PADDING 1
#define SIGMA 5.5
#define SHARPEN_AMOUNT 3
#define NOISE_THRESHOLD 2000
#define PI 3.14159265

//TODO cleaning

int cmpfunc (const void * a, const void * b) {
    return ( *(unsigned short*)a - *(unsigned short*)b );
}
void swap(short* a, short* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}
int partition (short arr[], short low, short high)
{
    short pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high- 1; j++)
    {

        if (arr[j] < pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

int bswap(int x)
{
    int mostleft;
    int mostright;

    mostleft=(x & 0x00ff) << 8;
    mostright=(x & 0xff00) >> 8;

    int result= mostleft |  mostright;
    return result;

}
void GaussFilterCreation(int kernel[][KERN_SIZE])
{
    // intialising standard deviation to 1.0

    double r, s = 2.0 * SIGMA * SIGMA;
    double GKernel[KERN_SIZE][KERN_SIZE];
    // sum is for normalization
    double sum = 0.0;

    // generating 5x5 kernel
    for (int x = -KERN_SIZE/2; x <= KERN_SIZE/2; x++) {
        for (int y = -KERN_SIZE/2; y <= KERN_SIZE/2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + KERN_SIZE/2][y + KERN_SIZE/2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + KERN_SIZE/2][y + KERN_SIZE/2];
        }
    }
    for (int i = 0; i < KERN_SIZE; ++i)
        for (int j = 0; j < KERN_SIZE; ++j)
            GKernel[i][j] /= sum;

    for (int i = 0; i < KERN_SIZE; ++i) {
        for (int j = 0; j < KERN_SIZE; ++j) {
            GKernel[i][j] =GKernel[i][j]*pow(2,(KERN_SIZE/2)+4);
            kernel[i][j]=(int)(GKernel[i][j]+0.5);
           // kernel[i][j]=GKernel[i][j];
        }
    }
}

void GaussFilterCreation2(double kernel[][KERN_SIZE])
{
    // intialising standard deviation to 1.0

    double r, s = 2.0 * SIGMA * SIGMA;
    double GKernel[KERN_SIZE][KERN_SIZE];
    // sum is for normalization
    double sum = 0.0;

    // generating 5x5 kernel
    for (int x = -KERN_SIZE/2; x <= KERN_SIZE/2; x++) {
        for (int y = -KERN_SIZE/2; y <= KERN_SIZE/2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + KERN_SIZE/2][y + KERN_SIZE/2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + KERN_SIZE/2][y + KERN_SIZE/2];
        }
    }
    for (int i = 0; i < KERN_SIZE; ++i)
        for (int j = 0; j < KERN_SIZE; ++j)
            GKernel[i][j] /= sum;

    for (int i = 0; i < KERN_SIZE; ++i) {
        for (int j = 0; j < KERN_SIZE; ++j) {
           // GKernel[i][j] =GKernel[i][j]*pow(2,(KERN_SIZE/2)+4);
            //kernel[i][j]=(int)(GKernel[i][j]+0.5);
             kernel[i][j]=GKernel[i][j];
        }
    }

}

int checksum(int kernel[][KERN_SIZE])
{
    int sum=0;
    for (int i = 0; i < KERN_SIZE; ++i)
        for (int j = 0; j < KERN_SIZE; ++j)
            sum+=kernel[i][j];
    return sum;
}
int main() {

    int bpp=2;
    int SIZE=WIDTH*HEIGHT*bpp;
    unsigned short pixel=0x0000;
    unsigned short blur_pixel=0x0000;
    int Gx[KERN_SIZE][KERN_SIZE]={1,0,-1,
                                  2,0,-2,
                                  1,0,-1};
    int Gy[KERN_SIZE][KERN_SIZE]={1,2,1,
                                  0,0,0,
                                  -1,-2,-1};
    int rows=0;
    FILE * pFile;
    long lSize,outlSize=(WIDTH+2*PADDING)*(HEIGHT+2*PADDING)*B;
    long rgbSize=WIDTH*HEIGHT*3;
    char * buffer;
    char * out_buffer;
    char * blur_buffer;
    char * sharp_buffer;
    char * edge_buffer;
    char * rgb_buffer;
    size_t result;
    size_t i=0,k=0;
    int j;
    int kernel[KERN_SIZE][KERN_SIZE];

    //IMG_1280x720_16bpp
    pFile = fopen ( "D:\\DemoApps\\intel5\\IMG_1280x720_16bpp.raw" , "rb" );
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}


    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    //lSize=SIZE;
    rewind (pFile);

    // allocate memory to contain the whole file:
    buffer = (char*) malloc (sizeof(char)*lSize);
    blur_buffer = (char*) malloc (sizeof(char)*lSize);
    sharp_buffer = (char*) malloc (sizeof(char)*lSize);
    edge_buffer = (char*) malloc (sizeof(char)*lSize);
    out_buffer = (char*) malloc (sizeof(char)*outlSize);
    rgb_buffer = (char*) malloc (sizeof(char)*rgbSize);

    if (buffer == NULL) {fputs ("Memory error reading buffer",stderr); exit (2);}
    if (out_buffer == NULL) {fputs ("Memory error padding buffer",stderr); exit (2);}
    if (blur_buffer == NULL) {fputs ("Memory error blur buffer",stderr); exit (2);}
    if (edge_buffer == NULL) {fputs ("Memory error edge buffer",stderr); exit (2);}
    if (sharp_buffer == NULL) {fputs ("Memory error sharp buffer",stderr); exit (2);}

    //Y*(width*B)+(B*X)

    result = fread (buffer,1,lSize,pFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}


    //PADDING
    while(rows<HEIGHT+2*PADDING) {

        if (rows < PADDING)
        {

            //left top corner
            for (i = i, j = 0; j < PADDING; i += 2, j++) {
                out_buffer[i] = buffer[0];
                out_buffer[i + 1] = buffer[1];
            }

            //above
            memcpy(&out_buffer[i], &buffer[0], WIDTH * B);

            i = i + WIDTH * 2;
            k = k + WIDTH * 2;

            //right top corner
            for (i = i, j = 0; j < PADDING; i += 2, j++) {
                out_buffer[i] = buffer[0 * (WIDTH * B) + (B * WIDTH-1)];
                out_buffer[i + 1] = buffer[0 * (WIDTH * B) + (B * WIDTH-1) + 1];
            }
        }
        else
        {
            if (rows > HEIGHT-1+PADDING) {

                //left down corner
                for (i = i, j = 0; j < PADDING; i += 2, j++) {
                    out_buffer[i] = buffer[(HEIGHT-1) * (WIDTH * B) + (B * 0)];
                    out_buffer[i + 1] = buffer[(HEIGHT-1) * (WIDTH * B) + (B * 0) + 1];
                }

                //under
                memcpy(&out_buffer[i], &buffer[(HEIGHT-1) * (WIDTH * B) + (B * 0)], WIDTH * B);

                i = i + WIDTH * 2;
                k = k + WIDTH * 2;

                //right down corner
                for (i = i, j = 0; j < PADDING; i += 2, j++) {

                    out_buffer[i] = buffer[(HEIGHT-1) * (WIDTH * B) + (B * (WIDTH-1))];
                    out_buffer[i + 1] = buffer[(HEIGHT-1) * (WIDTH * B) + (B * (WIDTH-1)) + 1];

                }
            } else{

                //left
                for (i=i,j=0;j<PADDING;i+=2,j++) {
                    out_buffer[i] = buffer[(rows-PADDING)*(WIDTH*B)+(B*0)];
                    out_buffer[i+1]=buffer[(rows-PADDING)*(WIDTH*B)+(B*0)+1];
                }

                //middle
                memcpy(&out_buffer[i], &buffer[(rows-PADDING)*(WIDTH*B)+(B*0)], WIDTH*B);
                i = i + WIDTH * 2;
                k = k + WIDTH * 2;

                //right
                for(i=i,j=0;j<PADDING;i+=2,j++) {
                    out_buffer[i] = buffer[(rows-PADDING)*(WIDTH*B)+(B*(WIDTH-1))];
                    out_buffer[i+1]=buffer[(rows-PADDING)*(WIDTH*B)+(B*(WIDTH-1))+1];
                }
            }
        }
        rows++;
    }

    rows=0;
    i=0;
    int a,b;
    int sum=0x0000;
    int STRIDE=2*PADDING+WIDTH;

   // double convolution;

    //Y*(width*B)+(B*X)
    rows=0;//x=j   y=rows

    int xpad,ypad;
    unsigned short v[KERN_SIZE*KERN_SIZE];
    for(a=0;a<KERN_SIZE*KERN_SIZE;a++)
        v[a]=0x0000;
    int count=0;
    double dkernel[KERN_SIZE][KERN_SIZE];
    GaussFilterCreation(kernel);
    GaussFilterCreation2(dkernel);
    sum=checksum(kernel);
    size_t contor=0;
    for (size_t i = 0,j=0; (i + 1) < lSize; i += 2,j++) {

        blur_buffer[i] = buffer[i];
        blur_buffer[i + 1] = buffer[i + 1];

        if (i > 2 * WIDTH * (rows + 1)) {
            rows++;
            j = 0;
        }
        xpad = j + PADDING;
        ypad = rows + PADDING;
        //pixel_coordinate_new = ((PADDING + rows - 1) * STRIDE * B + (PADDING + j) * B);
        count = 0;

        //j==418&&rows==383
        // if(j==418&&rows==383) {
             //printf("%d \n",sum);
             int convolution=0;
             int convolutionX=0;
             int convolutionY=0;
              //printf("%d  %d \n ",xpad ,xpad + KERN_SIZE / 2);
              for (int n = ypad - KERN_SIZE / 2, b = 0; n <= ypad + KERN_SIZE / 2; n++, b++) {
                  for (int m = xpad - KERN_SIZE / 2, a = 0; m <= xpad + KERN_SIZE / 2; m++, a++) {
                      int coord = ((n) * STRIDE * B + (m) * B);
                      //printf("%x%x ", out_buffer[coord] & 0xff, out_buffer[coord + 1] & 0xff);
                      //kernel[a][b]= (out_buffer[coord] << 8) + (out_buffer[coord + 1]);
                      //kernel[a][b]=bswap(kernel[a][b]);

                      pixel =(out_buffer[coord] << 8) + (out_buffer[coord + 1]);
                      pixel=bswap(pixel);

                      //gaussian blur:
                      //convolution=convolution+ pixel* (short)kernel[a][b];
                      //gaussian 2:
                      convolution=convolution+ pixel*dkernel[a][b];

                      //median filter:
                     // v[count] = (out_buffer[coord] << 8) + (out_buffer[coord + 1]);
                     // v[count] = bswap(v[count]);

                     //EDGE DETECTION

                      convolutionX=convolutionX+pixel*Gx[a][b];
                      convolutionY=convolutionY+pixel*Gy[a][b];

                      count++;

                  }
                // printf("\n");

              }

              //median filter
             // qsort(v, KERN_SIZE * KERN_SIZE, B, cmpfunc);
             // pixel = v[KERN_SIZE * KERN_SIZE / 2];


             //GAUSS:
             pixel=(short)convolution;

             blur_buffer[i] = pixel & 0xFF;
             blur_buffer[i + 1] = (pixel & 0xFF00) >> 8;




            //edge detection
            int test=(int)sqrt(convolutionX*convolutionX+convolutionY*convolutionY);
            unsigned short G=(short)sqrt(convolutionX*convolutionX+convolutionY*convolutionY);
                pixel = G&0xffff;



                //REMEMBER KERNEL 3 and PADDING 1 FOR THIS OR IT DOESN'T WORK
            //normal edge
            edge_buffer[i] = pixel & 0xFF;
            edge_buffer[i + 1] = (pixel & 0xFF00) >> 8;

            //RGB edges
            if(convolutionX==0)
                convolutionX=1;
            double angle = atan2(convolutionY, convolutionX);
            int res = (angle * 180) / PI;
            if (res < 0)
                res = res+360;

           // printf("%d\n",res);
            unsigned char red,green,blue;

            if(G>6000) {
                G=G/129;
               // if(G>=255)
                 //   G=255;
                G=G&0xff;
                if(res>0&&res<=30)
                {
                    red=G&0xff;
                    green=0;
                    blue=0;
                }
                if(res>30&&res<=60)
                {
                    red=G&0xff;
                    green=G&0xff;
                    blue=0;
                }
                if(res>60&&res<=90)
                {
                    red=0;
                    green=G&0xff;
                    blue=0;
                }
                if(res>90&&res<=120)
                {
                    red=0;
                    green=G&0xff;
                    blue=G&0xff;
                }
                if(res>120&&res<=150)
                {
                    red=0;
                    green=0;
                    blue=G&0xff;
                }
                if(res>150&&res<=180)
                {
                    red=G&0xff;
                    green=0;
                    blue=G&0xff;
                }
                if(res>180&&res<=210)
                {
                    red=G&0xff;
                    green=0;
                    blue=0;
                }
                if(res>210&&res<=240)
                {
                    red=G&0xff;
                    green=G&0xff;
                    blue=0;
                }
                if(res>240&&res<=270)
                {
                    red=0;
                    green=G&0xff;
                    blue=0;
                }
                if(res>270&&res<=300)
                {
                    red=0;
                    green=G&0xff;
                    blue=G&0xff;
                }if(res>300&&res<=330)
                {
                    red=0;
                    green=0;
                    blue=G&0xff;
                }if(res>330&&res<=360)
                {
                    red=G&0xff;
                    green=0;
                    blue=G&0xff;
                }

                rgb_buffer[i+contor] = red;
                rgb_buffer[i+contor+1] =green;
                rgb_buffer[i+contor+2] = blue;
            }
           // {
            //    rgb_buffer[i+contor] = 0x00;
           //     rgb_buffer[i+contor+1] = 0x00;
           //     rgb_buffer[i+contor+2] = 0x00;
         //   }

            contor++;


         // }
    }


    unsigned short sharp_pixel;
    for (size_t i = 0; (i + 1) < lSize; i += 2)
    {

        //original image pixel
        pixel=(buffer[i] << 8) + (buffer[i + 1]);
        pixel=bswap(pixel);

        //blur_pixel=pixe
        //blurred image pixel
        blur_pixel=(blur_buffer[i] << 8) + (blur_buffer[i + 1]);
        blur_pixel=bswap(blur_pixel);

        //S= O + A(O-B), if |O-B| > T
        // S - sharpened image
        // O - original image
        // B - blurred image
        // A - sharpen amount
        // T - noise threshold

        if(abs(pixel-blur_pixel)>NOISE_THRESHOLD) {
            if ((pixel + SHARPEN_AMOUNT * (pixel - blur_pixel)) > 0xffff)
                pixel = 0xffff;
            else {
                if ((pixel + SHARPEN_AMOUNT * (pixel - blur_pixel)) < 0x0000)
                   pixel = 0x0000;
                else
                    pixel = (pixel + SHARPEN_AMOUNT * (pixel - blur_pixel));

            }

        }
       // pixel= SHARPEN_AMOUNT*(pixel-blur_pixel);
        //sharp_pixel=pixel-blur_pixel;

        sharp_buffer[i] = pixel & 0xFF;
        sharp_buffer[i + 1] = (pixel & 0xFF00) >> 8;
    }



    // Open File(s) for writing and write buffers
    FILE *f_dst = fopen("blur_image.raw", "wb");
    if(f_dst == NULL)
    {
       printf("ERROR - Failed to open file for writing\n");
        exit(1);
    }
    if(fwrite(blur_buffer, 1, lSize, f_dst) != lSize)
    {
        printf("ERROR - Failed to write %i bytes to file\n", lSize);
        exit(1);
    }

    FILE *f_dst2 = fopen("sharp_image.raw", "wb");
    if(f_dst2 == NULL)
    {
        printf("ERROR - Failed to open file for writing\n");
        exit(1);
    }
    if(fwrite(sharp_buffer, 1, lSize, f_dst2) != lSize)
    {
        printf("ERROR - Failed to write %i bytes to file\n", lSize);
        exit(1);
    }

    FILE *f_dst3 = fopen("edges_image.raw", "wb");
    if(f_dst3 == NULL)
    {
        printf("ERROR - Failed to open file for writing\n");
        exit(1);
    }
    if(fwrite(edge_buffer, 1, lSize, f_dst3) != lSize)
    {
        printf("ERROR - Failed to write %i bytes to file\n", lSize);
        exit(1);
   }

    FILE *f_dst4 = fopen("rgb_edges_image.raw", "wb");
    if(f_dst4 == NULL)
    {
        printf("ERROR - Failed to open file for writing\n");
        exit(1);
    }


    if(fwrite(rgb_buffer, 1, rgbSize, f_dst4) != rgbSize)
    {
        printf("ERROR - Failed to write %i bytes to file\n", rgbSize);
        exit(1);
    }


// Close Files
    fclose(f_dst);
    fclose(f_dst2);
    fclose(f_dst3);
    fclose(f_dst4);
    f_dst = NULL;
    f_dst2 = NULL;
    f_dst3 = NULL;
    f_dst4 = NULL;

    fclose (pFile);
    free(buffer);
    free(out_buffer);
    free(blur_buffer);
    free(sharp_buffer);
    free(rgb_buffer);

    return 0;
}
