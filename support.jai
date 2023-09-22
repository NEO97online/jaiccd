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

ccd_support_t :: struct {
    v: Vector3;  //!< Support point in minkowski sum
    v1: Vector3; //!< Support point in obj1
    v2: Vector3; //!< Support point in obj2
};

/**** INLINES ****/
ccdSupportCopy :: inline (d: *ccd_support_t, s: *ccd_support_t) {
    d.* = s.*;
}

/**
 * Computes support point of obj1 and obj2 in direction dir.
 * Support point is returned via supp.
 */
__ccdSupport :: (obj1: *void, obj2: *void, _dir: *Vector3, ccd: *ccd_t, supp: *ccd_support_t) {
    dir: Vector3;

    ccdVec3Copy(*dir, _dir);

    ccd.support1(obj1, *dir, *supp.v1);

    ccdVec3Scale(*dir, -CCD_ONE);
    ccd.support2(obj2, *dir, *supp.v2);

    ccdVec3Sub2(*supp.v, *supp.v1, *supp.v2);
}