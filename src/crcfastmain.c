/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "crcfast.h"


/**
 * Main
 */

int main(int argc, char* argv[]) __attribute__((weak));
int main(int argc, char* argv[]){
	clock_t   ts, te;
	double    td;
	uint32_t  ck;
	uint64_t  generator = 0x82f63b78U;
	CRCFAST   crc32posix;
	CRCFAST   crc32c;
	CRCFAST   crc32;
	size_t    bufLen;
	char*     buf;
	long long i;
	
	
	/* El-cheapo argument parsing */
	if(argc > 1){
		for(i=1;i<argc-1;i++){
			if(strcmp(argv[i], "--generator") == 0||
			   strcmp(argv[i], "-g"         ) == 0){
				if      (strcmp(argv[i+1], "posix"     ) == 0){
					generator = 0xedb88320U;
				}else if(strcmp(argv[i+1], "castagnoli") == 0){
					generator = 0x82f63b78U;
				}else{
					generator = strtoull(argv[i+1], 0, 0);
				}
				i++;
			}
		}
	}
	
	
	/* Printout a message to ensure we got the generator right. */
	printf("Using generator polynomial %08llx",
	       (unsigned long long)generator);
	if      (generator == 0xedb88320U){
		printf(" (posix)");
	}else if(generator == 0x82f63b78U){
		printf(" (castagnoli)");
	}
	printf(".\n");
	
	
	/* Create a bunch of CRC objects. */
	crcfastInit(&crc32posix, 32, 0xedb88320U);
	crcfastInit(&crc32c,     32, 0x82f63b78U);
	crcfastInit(&crc32,      32, generator);
	
	
	/**
	 * Create a test buffer to compute the CRC over, making it at least
	 * *somewhat* challenging (not all-zeros).
	 */
	bufLen = 1*1024*1024*1024;
	buf    = calloc(1, bufLen);
	for(i=0;(size_t)i<bufLen;i+=4096){
		buf[i] = 0;/* "Touch" the pages to ensure they're in. */
	}
	for(i=0;i<256;i++){
		buf[2039+i] = i;
	}
	
	
	/* Do the actual benchmarking and print the results. */
	ts = clock();
	ck = ~crcfast(&crc32, ~0, buf, bufLen);
	te = clock();
	td = (te-ts)*1./CLOCKS_PER_SEC;
	printf("%08x (time: %f, BW: %.3f GiB/s)\n",
	       ck, td, bufLen/td/(1024.*1024.*1024.));
	
	
	/* Cleanup */
	crcfastDelete(&crc32posix);
	crcfastDelete(&crc32c);
	crcfastDelete(&crc32);
	
	
	/* Return */
	return 0;
}
