CCD_INIT :: (ccd: *ccd_t) {
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

size_t :: u64;

#load "ccd.jai";