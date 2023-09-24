CCD_INIT :: (ccd: *CCD_Config) {
    ccd.first_dir = ccdFirstDirDefault;
    ccd.support1 = null;
    ccd.support2 = null;
    ccd.center1 = null;
    ccd.center2 = null;

    ccd.max_iterations = #run cast(u32)(0) - 1;
    ccd.epa_tolerance = 0.0001;
    ccd.mpr_tolerance = 0.0001;
    ccd.dist_tolerance = 0.000001;
}

#scope_module

#import "Basic";
#import "Sort";

size_t :: u64;

_ccd_prefetch :: (a: Any) {}

#load "ccd.jai";
#load "alloc.jai";
#load "list.jai";
#load "vec3.jai";
#load "simplex.jai";
#load "support.jai";
#load "polytope.jai";
#load "mpr.jai";