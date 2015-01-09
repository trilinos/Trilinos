/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
/**
 * @file
 * @author  H. Carter Edwards
 * @date    October 2002
 *
 * @par Product Registry
 *   These 'C' functions support a registry of products
 *   that are linked into an application.  Each product
 *   has a name, version, and optional qualifier.
 *   Re-registration of products is permitted; however,
 *   the second and subsequent registrations must be
 *   identical to the first registration.  The product
 *   registry may be searched and iterated.
 *
 * @par Why 'C' instead of 'C++'
 *   As a general utility this simple product registry may
 *   be called from C++, C, or even FORTRAN.
 */

#ifndef STK_UTIL_ENVIRONMENT_PRODUCT_REGISTRY_H
#define STK_UTIL_ENVIRONMENT_PRODUCT_REGISTRY_H
#include <stddef.h>

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * @brief Extern "C" function <b>product_registry_add</b> provides a means for c
 * programs to register a product.  The arguments are passed on to
 * <b>ProductRegistry::addProduct()</b>/
 *
 * @param name			a <b>char</b> const pointer to the product's
 *				name.
 */
extern void product_registry_add(const char *name);

/**
 * @brief Extern "C" function <b>product_registry_add_tpl</b> provides a means for c
 * programs to register a tpl.  The arguments are passed on to
 * <b>ProductRegistry::addTpl()</b>/
 *
 * @param name			a <b>char</b> const pointer to the product's
 *				name.
 *
 * @param version              a <b>char</b> const pointer to the product's
 *                             version string.
 *
 * @param qualifier            a <b>char</b> const pointer to the product's
 *                             qualifier string.
 */
extern void product_registry_add_tpl(const char *name, const char *version, const char *qualifier);

/**
 * @brief Extern "C" function <b>product_registry_size</b> returns the number of
 * products in the registry.
 *
 * @return			an <b>int</b> value of the number of products in the
 *				registry.
 */
extern size_t product_registry_size();

#if defined(__cplusplus)
} /* extern "C" */
#endif


#endif /* STK_UTIL_ENVIRONMENT_PRODUCT_REGISTRY_H */
