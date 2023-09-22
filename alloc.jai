/***
 * libccd
 * ---------------------------------
 * Copyright (c)2010 Daniel Fiser <danfis@danfis.cz>
 *
 *
 *  This file is part of libccd.
 *
 *  Distributed under the OSI-approved BSD License (the "License");
 *  see accompanying file BDS-LICENSE for details or see
 *  <http://www.opensource.org/licenses/bsd-license.php>.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the License for more information.
 */

/**
 * Functions and macros required for memory allocation.
 */

/* Memory allocation: */
__CCD_ALLOC_MEMORY :: (type: Type, ptr_old: *void, size: size_t) {
    cast(*type)realloc(ptr_old, (size));
}

/** Allocate memory for one element of type.  */
CCD_ALLOC :: (type: Type) #expand {
    __CCD_ALLOC_MEMORY(type, null, size_of(type))
}

/** Allocate memory for array of elements of type type.  */
CCD_ALLOC_ARR :: (type: Type, num_elements: int) #expand {
    __CCD_ALLOC_MEMORY(type, null, size_of(type) * num_elements)
}

CCD_REALLOC_ARR :: (ptr: *void, type: Type, num_elements: int) #expand {
    __CCD_ALLOC_MEMORY(type, ptr, size_of(type) * num_elements)
}
