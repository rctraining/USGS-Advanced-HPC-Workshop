/*
 * Copyright (C) 2015  Timothy Brown
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ATTS_H
#define ATTS_H

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef _POSIX_C_SOURCE
#  undef _POSIX_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L
#define _GNU_SOURCE

#ifdef __BIGGEST_ALIGNMENT__
#  define ALIGNMENT __BIGGEST_ALIGNMENT__
#else
#  define ALIGNMENT 64
#endif

#if defined(__INTEL__)
#  define ATT_ALIGN __attribute__((aligned(ALIGNMENT)))
#  define ASSUME_ALIGNED __assume_aligned
#elif defined(__GNUC__)
#  define ATT_ALIGN __attribute__((aligned(ALIGNMENT)))
#  define ASSUME_ALIGNED __builtin_assume_aligned
#else
#  define ATT_ALIGN
#  define ASSUME_ALIGNED
#endif


#ifdef __cplusplus
}                               /* extern "C" */
#endif

#endif                          /* ATTS_H */
