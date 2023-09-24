# jaiccd
Jai port of https://github.com/danfis/libccd

## Features

- [x] `ccdGJKIntersect`
- [x] `ccdGJKPenetration`
	- Working as long as the EPA tolerance isn't too low. At `0.0001` tolerance, getting an infinite loop. Maybe this is because of 32-bit float precision?
	- Seems to be giving the wrong position. Not sure if this is also an upstream libccd issue...
	- If you need penetration, better to use `ccdMPRPenetration` for now. 
- [x] `ccdMPRIntersect`
- [x] `ccdMPRPenetration`

## Help Wanted

PRs are very welcome, especially for the following things:
- MPR implementations
- Porting of test suites