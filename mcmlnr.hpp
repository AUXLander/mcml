#pragma once

/***********************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void nrerror(const char error_text[]);