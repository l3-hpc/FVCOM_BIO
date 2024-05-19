# FVCOM_BIO

For these bio models, I am rewriting so all of them can use the same mod_bio_3D.F.

WARNING: this is work in progress.  In particular, since I am trying the above, when I update one model I may break another.  I will make a tagged version when everything is fine.

One more thing: Don't Use CGEM.  My code is in GitHub is so I don't lose work in progress, not because it's ready for prime time.

## FVCOM Licensing
The files mod_bio_3D.F, bio_mixing.F, mod_1D.F, mod_parameter.F, makefile, make.inc were originally from FVCOM.  There are some changes, but those are largely written by FVCOM folks.

Please read the [FVCOM Licencing agreement](https://github.com/FVCOM-GitHub/FVCOM/blob/main/LICENSE.md) before using anything.

## BIO_TP
BIO_TP is a simple total phosphorus model.  TP moves and sinks, with option to settle to the bottom or sink out.  It was originally written by taking Mark Rowe's GEM routine and carefully removing everything except the sinking.  'Carefully', because we wanted the same, so I was a bit overboard cautious, which made it a bit awkward.  This version cleans that up a bit.  Perhaps it will help with single/double precision errors...the sinking was calculated as (sink), but then converted to per unit time (sink/dt), and then multiplied by dt ('source*dt').  

## BIO_WQEM
A while back, I decided to try rewriting GOMDOM(now WQEM) so it can be called by FVCOM.  It compiled and ran, but that's as far as testing went.  It was just to show the idea of how to add WQEM, rather than actually finish it.  So now I've gone back to it, and trying to make it actually work.

### BTW: Remember the weirdness?
I was showing Wilson and Cody the FVCOM-GOMDOM code, and there was some really weird copying back and forth that didn't make sense.  I just said (shrug) 'it wasn't me'.  Now I realize the weirdness was because you actually needed to do that when you had phyto, zoo, det, etc. separated into categories(BIO_P, BIO_Z, etc.) and had to add them back to a main BIO_VAR.  My BIOs just have one 'group', so the copy was really doing nothing but wasting time. Anyway, this version cleans that up a bit. 

It is sort of like the above weird multiplying and dividing by dt in the BIO_TP model...you don't really notice unless you get rid of most of the code.  

## BIO_CGEM
Since CGEM was rewritten for SCHISM to be as separate as possible from the hydro, I think it can be easily used with FVCOM.  So I'm trying it out.  The current status is that it compiles and runs, but doesn't actually work.  To finish this, I'll probably look around for a smaller test grid.

## Citations
I'll put acknowledgements and citations here later.
