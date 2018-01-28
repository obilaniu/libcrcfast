/* Include Guards */
#ifndef INCLUDE_CRCFAST_H
#define INCLUDE_CRCFAST_H


/**
 * When building a library, it is a good idea to expose as few as possible
 * internal symbols (functions, objects, data structures). Not only does it
 * prevent users from relying on private portions of the library that are
 * subject to change without any notice, but it can have performance
 * advantages:
 * 
 *   - It can make shared libraries link faster at dynamic-load time.
 *   - It can make internal function calls faster by bypassing the PLT.
 * 
 * Thus, the compilation will by default hide all symbols, while the API
 * headers will explicitly mark public the few symbols the users are permitted
 * to use with a tag CRCFAST_PUBLIC. We also define a tag CRCFAST_HIDDEN, since
 * it may be required to explicitly tag certain C++ types as visible in order
 * for exceptions to function correctly.
 * 
 * Additional complexity comes from non-POSIX-compliant systems, which
 * artificially impose a requirement on knowing whether we are building or
 * using a DLL.
 * 
 * The above commentary and below code is inspired from
 *                   'https://gcc.gnu.org/wiki/Visibility'
 */

#if    defined(CRCFAST_SHAREDOBJECT) &&  defined(_WIN32) &&  defined(CRCFAST_BUILDING_DLL)
# define CRCFAST_PUBLIC __declspec(dllexport)
# define CRCFAST_HIDDEN
#elif  defined(CRCFAST_SHAREDOBJECT) &&  defined(_WIN32) && !defined(CRCFAST_BUILDING_DLL)
# define CRCFAST_PUBLIC __declspec(dllimport)
# define CRCFAST_HIDDEN
#elif  defined(CRCFAST_SHAREDOBJECT) &&  __GNUC__ >= 4
# define CRCFAST_PUBLIC __attribute__((visibility("default")))
# define CRCFAST_HIDDEN __attribute__((visibility("hidden")))
#else
# define CRCFAST_PUBLIC
# define CRCFAST_HIDDEN
#endif
# define CRCFAST_STATIC static



/* Includes */
#include <stdint.h>
#include <stddef.h>
#include <errno.h>



/* Defines */
#define CRCFAST_SUCCESS        (        0)
#define CRCFAST_EINVAL         (  -EINVAL)
#define CRCFAST_ENOMEM         (  -ENOMEM)



/**
 * Forward declarations & Typedefs
 */

typedef uint64_t CRCFAST_UINT;

struct CRCFAST;

typedef struct CRCFAST        CRCFAST;


/**
 * The CRCFAST Object.
 */

struct CRCFAST{
	const void* VTABLE;
	const void* LUT;
	uint64_t    generator;
	int         degree;
};


/* Function Prototypes */

/**
 * @brief Initialize new CRCFAST object.
 * 
 * @param [out] self      The initialized object
 * @param [in]  degree    The degree of the polynomial.
 *                        ***** WARNING: CURRENTLY ONLY SUPPORT degree==32 *****
 * @param [in]  generator The polynomial to be used, in bit-reversed form
 *                        without the leading power. Examples:
 * 
 *                            Degree |  Generator  | Name
 *                            -------+-------------+---------------------
 *                              32   | 0xedb88320U | CRC32  (Ethernet/POSIX/Zlib)
 *                              32   | 0x82f63b78U | CRC32C (Castagnoli)
 * 
 * @returns CRCFAST_SUCCESS if successful; A non-zero error code otherwise.
 */

CRCFAST_PUBLIC int          crcfastInit  (CRCFAST*       self,
                                          int            degree,
                                          uint64_t       generator);


/**
 * @brief Update a CRC remainder over a data buffer with the given CRC algorithm.
 * 
 * Does NOT handle either input or output inversions/reflections. This is the
 * responsibility of the caller.
 * 
 * @param [in] self    CRC algorithm object to use.
 * @param [in] crc     Initial CRC remainder
 * @param [in] buf     Pointer to data buffer
 * @param [in] bufLen  Length of data buffer, in bytes
 * @return The new CRC remainder.
 */

CRCFAST_PUBLIC CRCFAST_UINT crcfast      (const CRCFAST* self,
                                          CRCFAST_UINT   crc,
                                          const char*    buf,
                                          size_t         bufLen);


/**
 * @brief Delete the CRCFAST object.
 * 
 * Releases any dynamically-allocated resources the object may have been using.
 * After this function is called, the object must not be used.
 * 
 * If the user allocated the CRCFAST object itself on the heap, it is the
 * user's responsibility to free() it.
 * 
 * @param [in] self  The CRCFAST object to destroy.
 */

CRCFAST_PUBLIC void         crcfastDelete(CRCFAST*       self);



/* End Include Guards */
#endif

