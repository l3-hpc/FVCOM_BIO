***To use the last confirmed working version of BIO_TP, the one that was used for publishing [this paper](https://doi.org/10.1016/j.ecolmodel.2024.110718), get the first tagged version, [the release](https://github.com/l3-hpc/FVCOM_BIO/releases/tag/v1.0)***

# FVCOM_BIO

For these bio models, I am rewriting so all of them can use the *same* mod_bio_3D.F.

***WARNING***
- This is work in progress.  In particular, since I am trying the above, when I update one model I may break another.  I will make a tagged version when everything is fine. (My code is in GitHub is so I don't lose work in progress, not because it's ready for prime time.)
- On last modification, I broke the TP model.  So likely, nothing works.
- Why do I leave broken code in GitHub?   After this 'first pass', I believe it will be relatively straightforward for someone to take over the implementation (given time and funding), but it would be extremely difficult for someone (other than me) to implement WQEM+CGEM+TP with FVCOM from scratch.  Wanna take over?  Contact me and I'll walk you through it.
- I think I would not continue working on this without first identifying a smaller test case.  I'm testing using Lake Erie, and it is just too big(takes too long to run, large outputs) for efficient dev and test.


## FVCOM Licensing
The files mod_bio_3D.F, bio_mixing.F, mod_1D.F, mod_parameter.F, makefile, make.inc were originally from FVCOM.  There are some changes, but those are largely written by FVCOM folks.

Please read the [FVCOM Licencing agreement](https://github.com/FVCOM-GitHub/FVCOM/blob/main/LICENSE.md) before using anything.

## BIO_TP
BIO_TP is a simple total phosphorus model.  TP moves and sinks, with option to settle to the bottom or sink out.  It was originally written by taking Mark Rowe's GEM routine and removing everything except the sinking. This version is cleaned up a bit, now that I'm less worried about it being 'exactly the same'.

## BIO_WQEM
A while back, I decided to try rewriting GOMDOM(now WQEM) so it can be called by FVCOM.  It compiled and ran, but that's as far as testing went.  It was just to show the idea of how to add WQEM, rather than actually finish it.  So now I've gone back to it, and trying to make it actually work.

So far, it appears to work, meaning when I tried different things, the output did what I expected.
- PAR wasn't going to the bottom.  (or, like, in more than .5m of water).  I found the initial LOC was way high so changed that, and got more PAR.  Then I changed to KD=1 and got more PAR. Phytoplankton grows with more par.
- I changed the reference temperatures to be lower and got more 'activity'.  (Lake Erie in January is cold...)
- Sinking is more noticable in variables that sink faster.  I made it 'not sink out'.
Here's a super quick viz of [PAR and GRE](https://youtu.be/JEc_AALF7x4).

Sinking is done by the same subroutine as used in BIO_TP.


## BIO_CGEM
Since CGEM was rewritten for SCHISM to be as separate as possible from the hydro, I think it can be easily used with FVCOM.  So I'm trying it out.  The current status is that it compiles and runs, but doesn't actually work.  To finish this, I'll probably look around for a smaller test grid.  (Update: I think the problem was with precision, see below.  A new test run is waiting in the queue.  I really need a smaller test case...)

## Precision!
Use double_precision, because bio-geo vars don't move much in 10 seconds, and especially with CGEM units, you are adding a tiny df * dt.  Use single precision output, or VisIt will crash.  Just be happy and follow make.inc.

## Citations
I'll put acknowledgements and citations here later.
