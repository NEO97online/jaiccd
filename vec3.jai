#import "Math";

ccd_real_t :: float32;

CCD_REAL_MAX :: FLOAT32_INFINITY;

CCD_ONE  : float32 : 1.0;
CCD_ZERO : float32 : 0.0;

ccd_vec3_t :: Vector3;

CCD_EPS  :: 0.0000000001;
CCD_REAL :: (x: $T) -> ccd_real_t { return x; }
CCD_SQRT :: (x: $T) -> ccd_real_t { return sqrt(x); }
CCD_FABS :: (x: $T) -> ccd_real_t { return abs(x); }
CCD_FMIN :: min;
CCD_FMAX :: max;

ccd_vec3_origin := Vector3.{0, 0, 0};

ccdSign :: inline (val: ccd_real_t) -> int {
    if (ccdIsZero(val)){
        return 0;
    }else if (val < CCD_ZERO){
        return -1;
    }
    return 1;
}

ccdIsZero :: inline (val: ccd_real_t) -> bool {
    return CCD_FABS(val) < CCD_EPS;
}

ccdEq :: inline (_a: ccd_real_t, _b: ccd_real_t) -> bool {
    ab: ccd_real_t;
    a, b: ccd_real_t;

    ab = CCD_FABS(_a - _b);
    if (CCD_FABS(ab) < CCD_EPS)
        return true;

    a = CCD_FABS(_a);
    b = CCD_FABS(_b);
    if (b > a){
        return ab < CCD_EPS * b;
    }else{
        return ab < CCD_EPS * a;
    }
}


ccdVec3X :: inline (v: *Vector3) -> ccd_real_t {
    return v.x;
}

ccdVec3Y :: inline (v: *Vector3) -> ccd_real_t {
    return v.y;
}

ccdVec3Z :: inline (v: *Vector3) -> ccd_real_t {
    return v.z;
}

ccdVec3Eq :: inline (a: *Vector3, b: *Vector3) -> bool {
    return ccdEq(ccdVec3X(a), ccdVec3X(b))
            && ccdEq(ccdVec3Y(a), ccdVec3Y(b))
            && ccdEq(ccdVec3Z(a), ccdVec3Z(b));
}

ccdVec3Len2 :: inline (v: *Vector3) -> ccd_real_t {
    return ccdVec3Dot(v, v);
}

ccdVec3Dist2 :: inline (a: *Vector3, b: *Vector3) -> ccd_real_t {
    ab: Vector3;
    ccdVec3Sub2(*ab, a, b);
    return ccdVec3Len2(*ab);
}

ccdVec3Set :: inline (v: *Vector3, x: ccd_real_t, y: ccd_real_t, z: ccd_real_t) {
    v.x = x;
    v.y = y;
    v.z = z;
}

ccdVec3Copy :: inline (v: *Vector3, w: *Vector3) {
	v.* = w.*;
}

ccdVec3Sub :: inline (v: *Vector3, w: *Vector3) {
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
}
ccdVec3Sub2 :: inline (d: *Vector3, v: *Vector3, w: *Vector3) {
    d.x = v.x - w.x;
    d.y = v.y - w.y;
    d.z = v.z - w.z;
}

ccdVec3Add :: inline (v: *Vector3, w: *Vector3) {
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
}

ccdVec3Scale :: inline (d: *Vector3, k: ccd_real_t) {
	d.* *= k;
}

ccdVec3Normalize :: inline (d: *Vector3) {
	k := 1.0 / sqrt(ccdVec3Len2(d));
	ccdVec3Scale(d, k);
}

ccdVec3Dot :: inline (a: *Vector3, b: *Vector3) -> ccd_real_t {
    dot := a.x * b.x;
    dot += a.y * b.y;
    dot += a.z * b.z;
    return dot;
}

ccdVec3Cross :: inline (d: *Vector3, a: *Vector3, b: *Vector3) {
    d.x = (a.y * b.z) - (a.z * b.y);
    d.y = (a.z * b.x) - (a.x * b.z);
    d.z = (a.x * b.y) - (a.y * b.x);
}

__ccdVec3PointSegmentDist2 :: inline (P: *Vector3, x0: *Vector3, b: *Vector3, witness: *Vector3) -> ccd_real_t {
    // The computation comes from solving equation of segment:
    //      S(t) = x0 + t.d
    //          where - x0 is initial point of segment
    //                - d is direction of segment from x0 (|d| > 0)
    //                - t belongs to <0, 1> interval
    // 
    // Than, distance from a segment to some point P can be expressed:
    //      D(t) = |x0 + t.d - P|^2
    //          which is distance from any point on segment. Minimization
    //          of this function brings distance from P to segment.
    // Minimization of D(t) leads to simple quadratic equation that's
    // solving is straightforward.
    //
    // Bonus of this method is witness point for free.

    dist, t: ccd_real_t;
    d, a: Vector3;

    // direction of segment
    ccdVec3Sub2(*d, b, x0);

    // precompute vector from P to x0
    ccdVec3Sub2(*a, x0, P);

    t  = -CCD_REAL(1.) * ccdVec3Dot(*a, *d);
    t /= ccdVec3Len2(*d);

    if (t < CCD_ZERO || ccdIsZero(t)){
        dist = ccdVec3Dist2(x0, P);
        if (witness)
            ccdVec3Copy(witness, x0);
    }else if (t > CCD_ONE || ccdEq(t, CCD_ONE)){
        dist = ccdVec3Dist2(b, P);
        if (witness)
            ccdVec3Copy(witness, b);
    }else{
        if (witness){
            ccdVec3Copy(witness, *d);
            ccdVec3Scale(witness, t);
            ccdVec3Add(witness, x0);
            dist = ccdVec3Dist2(witness, P);
        }else{
            // recycling variables
            ccdVec3Scale(*d, t);
            ccdVec3Add(*d, *a);
            dist = ccdVec3Len2(*d);
        }
    }

    return dist;
}

ccdVec3PointSegmentDist2 :: (P: *Vector3, x0: *Vector3, b: *Vector3, witness: *Vector3) -> ccd_real_t {
    return __ccdVec3PointSegmentDist2(P, x0, b, witness);
}

ccdVec3PointTriDist2 :: (P: *Vector3, x0: *Vector3, B: *Vector3, C: *Vector3, witness: *Vector3) -> ccd_real_t {
    // Computation comes from analytic expression for triangle (x0, B, C)
    //      T(s, t) = x0 + s.d1 + t.d2, where d1 = B - x0 and d2 = C - x0 and
    // Then equation for distance is:
    //      D(s, t) = | T(s, t) - P |^2
    // This leads to minimization of quadratic function of two variables.
    // The solution from is taken only if s is between 0 and 1, t is
    // between 0 and 1 and t + s < 1, otherwise distance from segment is
    // computed.

    d1, d2, a: Vector3;
    u, v, w, p, q, r, d: ccd_real_t;
    s, t, dist, dist2: ccd_real_t;
    witness2: Vector3;

    ccdVec3Sub2(*d1, B, x0);
    ccdVec3Sub2(*d2, C, x0);
    ccdVec3Sub2(*a, x0, P);

    u = ccdVec3Dot(*a, *a);
    v = ccdVec3Dot(*d1, *d1);
    w = ccdVec3Dot(*d2, *d2);
    p = ccdVec3Dot(*a, *d1);
    q = ccdVec3Dot(*a, *d2);
    r = ccdVec3Dot(*d1, *d2);

    d = w * v - r * r;
    if (ccdIsZero(d)){
        // To avoid division by zero for zero (or near zero) area triangles
        s, t = -1.;
    }else{
        s = (q * r - w * p) / d;
        t = (-s * r - q) / w;
    }

    if ((ccdIsZero(s) || s > CCD_ZERO)
            && (ccdEq(s, CCD_ONE) || s < CCD_ONE)
            && (ccdIsZero(t) || t > CCD_ZERO)
            && (ccdEq(t, CCD_ONE) || t < CCD_ONE)
            && (ccdEq(t + s, CCD_ONE) || t + s < CCD_ONE)){

        if (witness){
            ccdVec3Scale(*d1, s);
            ccdVec3Scale(*d2, t);
            ccdVec3Copy(witness, x0);
            ccdVec3Add(witness, *d1);
            ccdVec3Add(witness, *d2);

            dist = ccdVec3Dist2(witness, P);
        }else{
            dist  = s * s * v;
            dist += t * t * w;
            dist += CCD_REAL(2.) * s * t * r;
            dist += CCD_REAL(2.) * s * p;
            dist += CCD_REAL(2.) * t * q;
            dist += u;
        }
    }else{
        dist = __ccdVec3PointSegmentDist2(P, x0, B, witness);

        dist2 = __ccdVec3PointSegmentDist2(P, x0, C, *witness2);
        if (dist2 < dist){
            dist = dist2;
            if (witness)
                ccdVec3Copy(witness, *witness2);
        }

        dist2 = __ccdVec3PointSegmentDist2(P, B, C, *witness2);
        if (dist2 < dist){
            dist = dist2;
            if (witness)
                ccdVec3Copy(witness, *witness2);
        }
    }

    return dist;
}