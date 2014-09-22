/*
 * Copyright (c) 2013, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

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

#ifndef STK_UTIL_REGISTRY_PRODUCT_REGISTRY_H
#define STK_UTIL_REGISTRY_PRODUCT_REGISTRY_H
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


#endif /* STK_UTIL_REGISTRY_PRODUCT_REGISTRY_H */
