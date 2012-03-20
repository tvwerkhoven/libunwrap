/*
 debugprint.h -- Print stuff if DEBUG 
 Copyright (C) 2011 Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef HAVE_DEBUGPRINT_H
#define HAVE_DEBUGPRINT_H

#ifdef HAVE_DEBUGPRINT

// Use like: DEBUGPRINT("%ju, %Lg\n", t_epoch.i, t_epoch.f);
#define DEBUGPRINT(fmt, ...) \
do { if (1) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
__LINE__, __func__, __VA_ARGS__); } while (0)


#else
#define DEBUGPRINT(fmt, ...) do { ; } while (0)
#endif



#endif // HAVE_DEBUGPRINT_H
