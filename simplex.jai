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

CCD_Simplex :: struct {
    ps: [4] CCD_Support;
    last: int; //!< index of last added point
}

ccdSimplexInit :: inline (s: *CCD_Simplex) {
    s.last = -1;
}

ccdSimplexSize :: inline (s: *CCD_Simplex) -> int {
    return s.last + 1;
}

ccdSimplexLast :: inline (s: *CCD_Simplex) -> *CCD_Support {
    return ccdSimplexPoint(s, s.last);
}

ccdSimplexPoint :: inline (s: *CCD_Simplex, idx: int) -> *CCD_Support {
    // here is no check on boundaries (JAI: what? this used a const pointer, i guess we need to do something here to nocheck boundary)
    return *s.ps[idx];
}
ccdSimplexPointW :: inline (s: *CCD_Simplex, idx: int) -> *CCD_Support {
    return *s.ps[idx];
}

ccdSimplexAdd :: inline (s: *CCD_Simplex, v: *CCD_Support) {
    // here is no check on boundaries in sake of speed
    s.last += 1;
    ccdSupportCopy(s.ps.data + s.last, v);
}

ccdSimplexSet :: inline (s: *CCD_Simplex, pos: size_t, a: *CCD_Support) {
    ccdSupportCopy(s.ps.data + pos, a);
}

ccdSimplexSetSize :: inline (s: *CCD_Simplex, size: int) {
    s.last = size - 1;
}

ccdSimplexSwap :: inline (s: *CCD_Simplex, pos1: size_t, pos2: size_t) {
    supp: CCD_Support;

    ccdSupportCopy(*supp, *s.ps[pos1]);
    ccdSupportCopy(*s.ps[pos1], *s.ps[pos2]);
    ccdSupportCopy(*s.ps[pos2], *supp);
}
