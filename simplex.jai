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

ccd_simplex_t :: struct {
    ps: [4] ccd_support_t;
    last: int; //!< index of last added point
}

ccdSimplexInit :: inline (s: *ccd_simplex_t) {
    s.last = -1;
}

ccdSimplexSize :: inline (s: *ccd_simplex_t) -> int {
    return s.last + 1;
}

ccdSimplexLast :: inline (s: *ccd_simplex_t) -> *ccd_support_t {
    return ccdSimplexPoint(s, s.last);
}

ccdSimplexPoint :: inline (s: *ccd_simplex_t, idx: int) -> *ccd_support_t {
    // here is no check on boundaries (JAI: what? this used a const pointer, i guess we need to do something here to nocheck boundary)
    return *s.ps[idx];
}
ccdSimplexPointW :: inline (s: *ccd_simplex_t, idx: int) -> *ccd_support_t {
    return *s.ps[idx];
}

ccdSimplexAdd :: inline (s: *ccd_simplex_t, v: *ccd_support_t) {
    // here is no check on boundaries in sake of speed
    s.last += 1;
    ccdSupportCopy(s.ps.data + s.last, v);
}

ccdSimplexSet :: inline (s: *ccd_simplex_t, pos: size_t, a: *ccd_support_t) {
    ccdSupportCopy(s.ps.data + pos, a);
}

ccdSimplexSetSize :: inline (s: *ccd_simplex_t, size: int) {
    s.last = size - 1;
}

ccdSimplexSwap :: inline (s: *ccd_simplex_t, pos1: size_t, pos2: size_t) {
    supp: ccd_support_t;

    ccdSupportCopy(*supp, *s.ps[pos1]);
    ccdSupportCopy(*s.ps[pos1], *s.ps[pos2]);
    ccdSupportCopy(*s.ps[pos2], *supp);
}
