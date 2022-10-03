#Dr. Ricardo Zamel's custom functions

#functions to convert from fold to log2ratio and vice versa FROM RICARDO
log2ratio2fold <- function (l2r) {
	options(warn = -1);
	fold <- ifelse(
		l2r >= 0,
		2^l2r,
		-1/2^l2r
	);
	return(fold);
}