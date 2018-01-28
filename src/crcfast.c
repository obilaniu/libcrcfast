/**
 * The CRC algorithm has the linearity property: That is, for
 * bitstrings A, B, Z=<zeroes> all of equal and arbitrary length, and
 * checksum bitstrings C, D, 0,
 * 
 * 
 *     CRC(C,   "" ) = C
 *     CRC(0,   Z  ) = 0
 *     CRC(0, A ^ B) = CRC(0, A) ^ CRC(0, B)
 *     CRC(C, A ^ B) = CRC(C, A) ^ CRC(0, B)
 *                   = CRC(C, B) ^ CRC(0, A)
 *     CRC(C, A ^ Z) = CRC(C, A) ^ CRC(0, Z)
 *                   = CRC(0, A) ^ CRC(C, Z)
 *                   = CRC(C, A)
 *     CRC(C, A || CRC(C, A)) = 0
 * 
 * The last one is true because the CRC of the message A with initialization C,
 * herein called D, is the remainder after long division of A * 2**32 by the
 * generator polynomial G; Therefore, if the bitstring -D == D is appended, the
 * new remainder is
 *     (-D+D)*(2**32) === 0*(2**32) === 0 (mod G)
 * 
 * 
 * This generalizes to wider divisions, which we exploit:
 * 
 *     CRC(init, abcdefghijklmnopabcdefghijklmnopabcdefghijklmnopabcdefghijklmnop)
 * =   CRC(init, abcd                                                            )
 *   ^ CRC(   0,     efgh                                                        )
 *   ^ CRC(   0,         ijkl                                                    )
 *   ^ CRC(   0,             mnop                                                )
 *   ^ CRC(   0,                 abcd                                            )
 *   ^ CRC(   0,                     efgh                                        )
 *   ^ CRC(   0,                         ijkl                                    )
 *   ^ CRC(   0,                             mnop                                )
 *   ^ CRC(   0,                                 abcd                            )
 *   ^ CRC(   0,                                     efgh                        )
 *   ^ CRC(   0,                                         ijkl                    )
 *   ^ CRC(   0,                                             mnop                )
 *   ^ CRC(   0,                                                 abcd            )
 *   ^ CRC(   0,                                                     efgh        )
 *   ^ CRC(   0,                                                         ijkl    )
 *   ^ CRC(   0,                                                             mnop)
 * 
 * 
 * If we track 16 independent CRC checksums, we can therefore work in 16-way
 * parallel fashion:
 * 
 *   //  Init
 *   crc0                             =init
 *   crc1=crc2=crc3=...=crcD=crcE=crcF=0
 *   //  Master Align Loop:
 *   while(p % 64 != 0):
 *       crc0 = CRC(crc0, *p++)
 *       bytesLeft--
 *   //  For large inputs:
 *   if(bytesLeft >= 128):
 *       // Prefetch
 *       blk0 = ((uint32_t*)p)[ 0]   // abcd
 *       blk1 = ((uint32_t*)p)[ 1]   // efgh
 *       blk2 = ((uint32_t*)p)[ 2]   // ijkl
 *       ...
 *       blkE = ((uint32_t*)p)[14]   // ijkl
 *       blkF = ((uint32_t*)p)[15]   // mnop
 *       //  Master Vector Loop:
 *       while(bytesLeft >= 128):
 *           //  XOR with 4 bytes:
 *           crc0 ^= blk0
 *           crc1 ^= blk1
 *           crc2 ^= blk2
 *           ...
 *           crcE ^= blk3
 *           crcF ^= blk3
 *           //  Pointer advance
 *           p         += 64
 *           bytesLeft -= 64
 *           //  Prefetch
 *           blk0 = ((uint32_t*)p)[ 0]   // abcd
 *           blk1 = ((uint32_t*)p)[ 1]   // efgh
 *           blk2 = ((uint32_t*)p)[ 2]   // ijkl
 *           ...
 *           blkE = ((uint32_t*)p)[14]   // ijkl
 *           blkF = ((uint32_t*)p)[15]   // mnop
 *           //  Skip CRC ahead by 64 bytes:
 *           crc0  = CRC(crc0, <64 zeros>)
 *           crc1  = CRC(crc1, <64 zeros>)
 *           crc2  = CRC(crc2, <64 zeros>)
 *           ...
 *           crcE  = CRC(crcE, <64 zeros>)
 *           crcF  = CRC(crcF, <64 zeros>)
 *       //  Master Vector Loop peeled exit:
 *       crc1 ^= CRC(crc0, ((uint32_t*)p)[ 0])
 *       crc2 ^= CRC(crc1, ((uint32_t*)p)[ 1])
 *       crc3 ^= CRC(crc2, ((uint32_t*)p)[ 2])
 *       ...
 *       crcE ^= CRC(crcD, ((uint32_t*)p)[13])
 *       crcF ^= CRC(crcE, ((uint32_t*)p)[14])
 *       crc0  = crcF
 *       p         += 60
 *       bytesLeft -= 60
 *   //  Master Dealign Loop:
 *   while(bytesLeft > 0):
 *       crc0 = CRC(crc0, *p++)
 *       bytesLeft--
 *   //  Return
 *   return crcF
 * 
 * However, we need a technique to do the skipahead by 64 zero bytes. This is, in
 * fact, doable efficiently; It corresponds to a multiplication by the polynomial
 * x**512 (mod G(x)), and can be implemented efficiently with PSHUFB, requiring
 * 4 byte-lookups for all 8 nibbles (4-bit groups) of the CRC-32 state plus 24 XORs.
 * 
 * 
 * PROBABLE LOOK AND TOTAL COST OF INNER LOOP:
 *   - 4 loads of 16 bytes on ports 2, 3
 *   - 4x16 -> 16x4 byte transpose:
 *        abcdefghijklmnop    aeimaeimaeimaeim
 *        abcdefghijklmnop -> bfjnbfjnbfjnbfjn
 *        abcdefghijklmnop -> cgkocgkocgkocgko
 *        abcdefghijklmnop    dhlpdhlpdhlpdhlp
 *     Note: ***Any permutation of the columns is acceptable, ESPECIALLY if it
 *              makes the transposition cheaper by being partial! The only
 *              important thing is for groups of 4 consecutive bytes in memory to
 *              be splayed across the same lanes of 4 XMM registers !!!!!!!! ***
 *     Probable cost: 12 shuffles on port 5
 *   - 4 PXORs with the current CRC state on ports 0, 1
 *   - 1 SHR & 1 AND on ports 0, 1
 *   - 32 aligned loads from LUTs, 32 PSHUFBs on port 5
 *   - 24 PXORs on ports 0, 1
 *   - Pointer increment on ports 1
 *   - Dec+Branch predicted taken on port 6
 * 
 * 
 * PROBABLE BOTTLENECK:
 *   - 32+12 shuffles on port 5
 * 
 * 
 * PROBABLE THROUGHPUT:
 *   - Estimated on Haswell+:
 *     44cc /  64 bytes (0.6875  cc/b; 1.45 b/cc)
 *   - If running on IvyBridge and below, **probably** double that to
 *     22cc /  64 bytes (0.34375 cc/b; 2.91 b/cc)
 *   - If using AVX2 VPSHUFB and hacking the above, **might** double that to
 *     44cc / 128 bytes (0.34375 cc/b; 2.91 b/cc)
 */

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "crcfast.h"
#include <immintrin.h>


/* Defines */
#define ALIGNED(align)   __attribute__((aligned(align)))
#define GETVTABLE(self)  ((const CRCFAST_VTABLE*)(self)->VTABLE)
#define XOR(a,b)         _mm_xor_si128 ((a), (b))
#define AND(a,b)         _mm_and_si128 ((a), (b))
#define SHIFTRIGHT(a,b)  _mm_srli_epi16((a), (b))
#define LOADA(p)         _mm_load_si128((const __m128i*)(p))



/* Data Structure Forward Declarations, Typedefs and Definitions */
struct CRCFAST_VTABLE;
struct CRC32_SSSE3_LUT;

typedef struct CRCFAST_VTABLE  CRCFAST_VTABLE;
typedef struct CRC32_SSSE3_LUT CRC32_SSSE3_LUT;



/**
 * @brief The virtual method table for CRC implementations.
 */

struct CRCFAST_VTABLE{
	CRCFAST_UINT (*update)(const CRCFAST* self,
	                       CRCFAST_UINT   crc,
	                       const char*    buf,
	                       size_t         bufLen);
	void         (*delete)(CRCFAST*       self);
};

/**
 * @brief The structure defining the SSSE3 CRC32 implementation's LUTs.
 */

struct CRC32_SSSE3_LUT{
	uint32_t S[256];
	uint8_t  K[32][16];
} ALIGNED(512);



/* Static Function Forward Declarations */
CRCFAST_STATIC CRCFAST_UINT crcfastUpdateCRC32Sarwate(const CRCFAST* self,
                                                      CRCFAST_UINT   crc,
                                                      const char*    buf,
                                                      size_t         bufLen);
CRCFAST_STATIC CRCFAST_UINT crcfastUpdateCRC32SSSE3  (const CRCFAST* self,
                                                      CRCFAST_UINT   crc,
                                                      const char*    buf,
                                                      size_t         bufLen);
CRCFAST_STATIC void         crcfastDeleteNoop        (CRCFAST*       self);
CRCFAST_STATIC void         crcfastDeleteDynamic     (CRCFAST*       self);
CRCFAST_STATIC void         crcfastInit32MakeSSSE3LUT(CRCFAST*       self);
CRCFAST_STATIC void         crcfastInit32DumpLUT     (const CRCFAST* self);
CRCFAST_STATIC int          crcfastInit32            (CRCFAST*       self,
                                                      uint32_t       generator);



/* Useful constants */
CRCFAST_STATIC const char           ZEROBYTES[64] = {0};
CRCFAST_STATIC const CRCFAST_VTABLE CRCFAST_CRC32_SSSE3_VTABLE = {
	.update = crcfastUpdateCRC32SSSE3,
	.delete = crcfastDeleteDynamic,
};
CRCFAST_STATIC const CRCFAST_VTABLE CRCFAST_CRC32_SSSE3_STATIC_VTABLE = {
	.update = crcfastUpdateCRC32SSSE3,
	.delete = crcfastDeleteNoop,
};
#include "crc32_generator82f63b78_ssse3_lut.h"
#include "crc32_generatoredb88320_ssse3_lut.h"


/* Static Function  */

/**
 * Update CRC32 on a fixed number of bytes using given CRCFAST object and
 * the Sarwate byte-at-a-time LUT.
 * 
 * Does not do any initial or final inversions.
 */

CRCFAST_STATIC CRCFAST_UINT crcfastUpdateCRC32Sarwate(const CRCFAST* self,
                                                      CRCFAST_UINT   crc,
                                                      const char*    buf,
                                                      size_t         bufLen){
	size_t         i;
	uint32_t       crc32 = crc;
	const uint8_t* ptr   = (const uint8_t*)buf;
	
	for(i=0;i<bufLen;i++){
		crc32 = (crc32>>8) ^ ((uint32_t*)self->LUT)[(crc32 ^ *ptr++) & 0xFFU];
	}
	
	return crc32;
}

/**
 * Update CRC32 on a fixed number of bytes using given CRCFAST object and
 * the SSSE3-accelerated codepath.
 * 
 * Does not do any initial or final inversions.
 */

CRCFAST_STATIC CRCFAST_UINT crcfastUpdateCRC32SSSE3  (const CRCFAST* self,
                                                      CRCFAST_UINT   crc,
                                                      const char*    buf,
                                                      size_t         bufLen){
	int misalign64 = (intptr_t)buf & 0x3F;
	
	if(bufLen < 128){
		return crcfastUpdateCRC32Sarwate(self, crc, buf, bufLen);
	}
	
	if(misalign64){
		misalign64 = 64-misalign64;
		crc        = crcfastUpdateCRC32Sarwate(self, crc, buf, misalign64);
		buf       += misalign64;
		bufLen    -= misalign64;
	}
	
	if(bufLen >= 128){
		const __m128i* const K = (const __m128i*)(&((CRC32_SSSE3_LUT*)self->LUT)->K);
		const __m128i  MASK0F  = _mm_set1_epi8(0x0F);
		const __m128i  TRANSP  = _mm_setr_epi8( 0,  4,  8, 12,
		                                        1,  5,  9, 13,
		                                        2,  6, 10, 14,
		                                        3,  7, 11, 15);
		
		uint32_t       crcv[16];
		__m128i        byte0, byte1, byte2, byte3;
		__m128i        blk0,  blk1,  blk2,  blk3;
		__m128i        tmp0,  tmp1,  tmp2,  tmp3;
		__m128i        nib0l, nib0h, nib1l, nib1h,
		               nib2l, nib2h, nib3l, nib3h;
		__m128i        xor00, xor01, xor02, xor03,
		               xor10, xor11, xor12, xor13,
		               xor20, xor21, xor22, xor23,
		               xor30, xor31, xor32, xor33,
		               xor40, xor41, xor42, xor43,
		               xor50, xor51, xor52, xor53,
		               xor60, xor61, xor62, xor63,
		               xor70, xor71, xor72, xor73;
		
		byte0 = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,crc >>  0);
		byte1 = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,crc >>  8);
		byte2 = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,crc >> 16);
		byte3 = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,crc >> 24);
		
		/* Prefetch next block */
		blk0  = LOADA(buf +  0);
		blk1  = LOADA(buf + 16);
		blk2  = LOADA(buf + 32);
		blk3  = LOADA(buf + 48);
		
		while(bufLen >= 128){
			blk0  = _mm_shuffle_epi8  (blk0,  TRANSP);
			blk1  = _mm_shuffle_epi8  (blk1,  TRANSP);
			blk2  = _mm_shuffle_epi8  (blk2,  TRANSP);
			blk3  = _mm_shuffle_epi8  (blk3,  TRANSP);
			tmp0  = _mm_unpacklo_epi32(blk0,  blk1);
			tmp1  = _mm_unpackhi_epi32(blk0,  blk1);
			tmp2  = _mm_unpacklo_epi32(blk2,  blk3);
			tmp3  = _mm_unpackhi_epi32(blk2,  blk3);
			blk0  = _mm_unpacklo_epi64(tmp0,  tmp2);
			blk1  = _mm_unpackhi_epi64(tmp0,  tmp2);
			blk2  = _mm_unpacklo_epi64(tmp1,  tmp3);
			blk3  = _mm_unpackhi_epi64(tmp1,  tmp3);
			
			byte0 = XOR               (byte0, blk0);
			byte1 = XOR               (byte1, blk1);
			byte2 = XOR               (byte2, blk2);
			byte3 = XOR               (byte3, blk3);
			
			buf    += 64;
			bufLen -= 64;
			
			/* Prefetch next block */
			blk0  = LOADA             (buf +  0);
			blk1  = LOADA             (buf + 16);
			blk2  = LOADA             (buf + 32);
			blk3  = LOADA             (buf + 48);
			
			nib0l = AND               (           byte0,     MASK0F);
			nib0h = AND               (SHIFTRIGHT(byte0, 4), MASK0F);
			nib1l = AND               (           byte1,     MASK0F);
			nib1h = AND               (SHIFTRIGHT(byte1, 4), MASK0F);
			nib2l = AND               (           byte2,     MASK0F);
			nib2h = AND               (SHIFTRIGHT(byte2, 4), MASK0F);
			nib3l = AND               (           byte3,     MASK0F);
			nib3h = AND               (SHIFTRIGHT(byte3, 4), MASK0F);
			xor00 = _mm_shuffle_epi8  (LOADA(&K[ 0]), nib0l);
			xor01 = _mm_shuffle_epi8  (LOADA(&K[ 1]), nib0l);
			xor02 = _mm_shuffle_epi8  (LOADA(&K[ 2]), nib0l);
			xor03 = _mm_shuffle_epi8  (LOADA(&K[ 3]), nib0l);
			xor10 = _mm_shuffle_epi8  (LOADA(&K[ 4]), nib0h);
			xor11 = _mm_shuffle_epi8  (LOADA(&K[ 5]), nib0h);
			xor12 = _mm_shuffle_epi8  (LOADA(&K[ 6]), nib0h);
			xor13 = _mm_shuffle_epi8  (LOADA(&K[ 7]), nib0h);
			xor20 = _mm_shuffle_epi8  (LOADA(&K[ 8]), nib1l);
			xor21 = _mm_shuffle_epi8  (LOADA(&K[ 9]), nib1l);
			xor22 = _mm_shuffle_epi8  (LOADA(&K[10]), nib1l);
			xor23 = _mm_shuffle_epi8  (LOADA(&K[11]), nib1l);
			xor30 = _mm_shuffle_epi8  (LOADA(&K[12]), nib1h);
			xor31 = _mm_shuffle_epi8  (LOADA(&K[13]), nib1h);
			xor32 = _mm_shuffle_epi8  (LOADA(&K[14]), nib1h);
			xor33 = _mm_shuffle_epi8  (LOADA(&K[15]), nib1h);
			xor40 = _mm_shuffle_epi8  (LOADA(&K[16]), nib2l);
			xor41 = _mm_shuffle_epi8  (LOADA(&K[17]), nib2l);
			xor42 = _mm_shuffle_epi8  (LOADA(&K[18]), nib2l);
			xor43 = _mm_shuffle_epi8  (LOADA(&K[19]), nib2l);
			xor50 = _mm_shuffle_epi8  (LOADA(&K[20]), nib2h);
			xor51 = _mm_shuffle_epi8  (LOADA(&K[21]), nib2h);
			xor52 = _mm_shuffle_epi8  (LOADA(&K[22]), nib2h);
			xor53 = _mm_shuffle_epi8  (LOADA(&K[23]), nib2h);
			xor60 = _mm_shuffle_epi8  (LOADA(&K[24]), nib3l);
			xor61 = _mm_shuffle_epi8  (LOADA(&K[25]), nib3l);
			xor62 = _mm_shuffle_epi8  (LOADA(&K[26]), nib3l);
			xor63 = _mm_shuffle_epi8  (LOADA(&K[27]), nib3l);
			xor70 = _mm_shuffle_epi8  (LOADA(&K[28]), nib3h);
			xor71 = _mm_shuffle_epi8  (LOADA(&K[29]), nib3h);
			xor72 = _mm_shuffle_epi8  (LOADA(&K[30]), nib3h);
			xor73 = _mm_shuffle_epi8  (LOADA(&K[31]), nib3h);
			byte0 = XOR(XOR(XOR(xor00, xor10), XOR(xor20, xor30)),
			            XOR(XOR(xor40, xor50), XOR(xor60, xor70)));
			byte1 = XOR(XOR(XOR(xor01, xor11), XOR(xor21, xor31)),
			            XOR(XOR(xor41, xor51), XOR(xor61, xor71)));
			byte2 = XOR(XOR(XOR(xor02, xor12), XOR(xor22, xor32)),
			            XOR(XOR(xor42, xor52), XOR(xor62, xor72)));
			byte3 = XOR(XOR(XOR(xor03, xor13), XOR(xor23, xor33)),
			            XOR(XOR(xor43, xor53), XOR(xor63, xor73)));
		}
		
		tmp0  = _mm_unpacklo_epi8 (byte0, byte1);
		tmp1  = _mm_unpackhi_epi8 (byte0, byte1);
		tmp2  = _mm_unpacklo_epi8 (byte2, byte3);
		tmp3  = _mm_unpackhi_epi8 (byte2, byte3);
		byte0 = _mm_unpacklo_epi16(tmp0,  tmp2);
		byte1 = _mm_unpackhi_epi16(tmp0,  tmp2);
		byte2 = _mm_unpacklo_epi16(tmp1,  tmp3);
		byte3 = _mm_unpackhi_epi16(tmp1,  tmp3);
		
		_mm_storeu_si128((__m128i*)&crcv[ 0], byte0);
		_mm_storeu_si128((__m128i*)&crcv[ 4], byte1);
		_mm_storeu_si128((__m128i*)&crcv[ 8], byte2);
		_mm_storeu_si128((__m128i*)&crcv[12], byte3);
		
		crcv[0x1] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x0], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x2] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x1], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x3] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x2], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x4] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x3], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x5] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x4], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x6] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x5], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x7] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x6], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x8] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x7], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0x9] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x8], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xA] ^= crcfastUpdateCRC32Sarwate(self, crcv[0x9], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xB] ^= crcfastUpdateCRC32Sarwate(self, crcv[0xA], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xC] ^= crcfastUpdateCRC32Sarwate(self, crcv[0xB], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xD] ^= crcfastUpdateCRC32Sarwate(self, crcv[0xC], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xE] ^= crcfastUpdateCRC32Sarwate(self, crcv[0xD], buf, 4); buf+= 4; bufLen -= 4;
		crcv[0xF] ^= crcfastUpdateCRC32Sarwate(self, crcv[0xE], buf, 4); buf+= 4; bufLen -= 4;
		crc        = crcv[0xF];
	}
	
	return crcfastUpdateCRC32Sarwate(self, crc, buf, bufLen);
}


/**
 * @brief A no-op destructor, for use by the objects with statically-allocated LUTs.
 */

CRCFAST_STATIC void         crcfastDeleteNoop        (CRCFAST*       self){
	(void)self;
}

/**
 * @brief A destructor for objects with dynamically-allocated LUTs
 */

CRCFAST_STATIC void         crcfastDeleteDynamic     (CRCFAST*       self){
	free((void*)self->LUT);
	
	/**
	 * We don't free the VTABLE since this software doesn't dynamically
	 * allocate VTABLEs.
	 */
}

/**
 * @brief Initialize SSSE3 CRC32 LUTs.
 * @param [in/out]  self  CRCFAST object whose CRC32 LUTs must be
 *                        initialized.
 */

CRCFAST_STATIC void         crcfastInit32MakeSSSE3LUT(CRCFAST*       self){
	CRC32_SSSE3_LUT* const LUT       = (CRC32_SSSE3_LUT*)self->LUT;
	uint32_t               generator = self->generator;
	int                    i, j;
	uint32_t               b;
	
	
	memset(LUT, 0, sizeof(*LUT));
	
	/**
	 * Compute S (Sarwate) Table for byte-at-a-time processing.
	 */
	
	for(b=0;b<256;b++){
		/**
		 * Walk from MSB to LSB of polynomial doing long-division by generator.
		 * 
		 * Remember that in BR ordering (which is what this code assumes),
		 * the MSB of the polynomial is the LSB of the byte and vice-versa.
		 */
		
		LUT->S[b] = b;
		for(i=0;i<8;i++){
			if(LUT->S[b] & 1){
				LUT->S[b] >>= 1;
				LUT->S[b]  ^= generator;
			}else{
				LUT->S[b] >>= 1;
			}
		}
	}
	
	
	/**
	 * Compute K (SSSE3) Table for 64-zero-byte sKipahead using S Table.
	 * 
	 * For each of the 16 possible values of each of the 8 nibbles in a CRC-32,
	 * compute the 64-zero-byte skip-ahead value, and splay its bytes 4-way.
	 */
	
	for(i=0;i<8;i++){
		for(j=0;j<16;j++){
			b = crcfastUpdateCRC32Sarwate(self, (uint32_t)j << i*4, ZEROBYTES, 64);
			LUT->K[4*i+0][j] = b >>  0;
			LUT->K[4*i+1][j] = b >>  8;
			LUT->K[4*i+2][j] = b >> 16;
			LUT->K[4*i+3][j] = b >> 24;
		}
	}
	
	/* Silences unused warning for this function. */
	(void)crcfastInit32DumpLUT;
}

/**
 * @brief Dump LUTs for this CRC32 object.
 * 
 * THIS IS INTERNAL UTILITY CODE.
 * 
 * Dump the LUT's contents as an #include'able C header.
 * 
 * @param [in]  self  The object for which the LUTs should be dumped.
 */

CRCFAST_STATIC void         crcfastInit32DumpLUT     (const CRCFAST* self){
	const CRC32_SSSE3_LUT* const LUT = (const CRC32_SSSE3_LUT*)self->LUT;
	const uint32_t         generator = self->generator;
	int                    i, j;
	
	char headerPath[FILENAME_MAX];
	sprintf(headerPath, "crc%d_generator%08llx_ssse3_lut.h",
	        self->degree,
	        (unsigned long long)generator);
	
	FILE* fp = fopen(headerPath, "w");
	if(!fp){
		fprintf(stderr, "Error: Failed to open %s\n", headerPath);
		exit(-1);
	}
	
	fprintf(fp, "const struct CRC%d_SSSE3_LUT crc%d_generator%08llx_ssse3_lut = {\n",
	            self->degree, self->degree, (unsigned long long)generator);
	fprintf(fp, "{\n\t");
	for(i=0;i<256;i++){
		const char* term = ", ";
		if(i == 255){
			term = ",\n";
		}else if((i&7) == 7){
			term = ",\n\t";
		}
		fprintf(fp, "0x%08x%s", LUT->S[i], term);
	}
	fprintf(fp, "},\n");
	fprintf(fp, "{\n");
	for(i=0;i<32;i++){
		fprintf(fp, "\t{");
		for(j=0;j<16;j++){
			fprintf(fp, "0x%02x%s", (uint8_t)LUT->K[i][j],
			        j==15 ? "" : ", ");
		}
		fprintf(fp, "},\n");
	}
	fprintf(fp, "}\n");
	fprintf(fp, "};\n");
	fflush(fp);
	fclose(fp);
}

/**
 * @brief Initialize a CRC32 checker from its generator polynomial.
 * @param [out] self       The initialized object.
 * @param [in]  generator  The generator polynomial to be used.
 * @return 
 */

CRCFAST_STATIC int          crcfastInit32            (CRCFAST*       self,
                                                      uint32_t       generator){
	/* If we have precomputed that LUT, give it out. */
	switch(generator){
		#define G(g)                                                  \
		    case 0x##g:                                               \
		        self->LUT    = &crc32_generator##g##_ssse3_lut;       \
		        self->VTABLE = &CRCFAST_CRC32_SSSE3_STATIC_VTABLE;    \
		        return CRCFAST_SUCCESS;
		
		G(edb88320)
		G(82f63b78)
		
		#undef G
	}
	
	/* Otherwise, dynamically allocate and initialize the LUT. */
	self->VTABLE = &CRCFAST_CRC32_SSSE3_VTABLE;
	self->LUT    = calloc(1, sizeof(CRC32_SSSE3_LUT));
	if(!self->LUT){
		return CRCFAST_ENOMEM;
	}
	crcfastInit32MakeSSSE3LUT(self);
	
	return CRCFAST_SUCCESS;
}


/**
 * Public Function Definitions
 */

CRCFAST_PUBLIC int          crcfastInit  (CRCFAST*       self,
                                          int            degree,
                                          uint64_t       generator){
	memset(self, 0, sizeof(*self));
	
	self->degree    = degree;
	self->generator = generator;
	
	if(degree == 32){
		return crcfastInit32(self, generator);
	}else{
		return CRCFAST_EINVAL;
	}
}
CRCFAST_PUBLIC CRCFAST_UINT crcfast      (const CRCFAST* self,
                                          CRCFAST_UINT   crc,
                                          const char*    buf,
                                          size_t         bufLen){
	return GETVTABLE(self)->update(self, crc, buf, bufLen);
}
CRCFAST_PUBLIC void         crcfastDelete(CRCFAST*       self){
	GETVTABLE(self)->delete(self);
	memset(self, 0, sizeof(*self));
}

