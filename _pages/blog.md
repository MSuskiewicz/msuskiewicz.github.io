---
layout: archive
title: "Blog"
permalink: /blog/
author_profile: true
redirect_from:
  - /blog
---

{% include base_path %}

Chances are that this blog will die at some point. It will be regularly updated initially, and less and less as the time goes by. Alas, that's what usually happens. But it's a risk I'm willing to take. I will write here short notes, mostly on scientific articles I've read.

Protein homodimerisation (20/09/2023)
------
A lot of proteins we see in nature are homodimers in which the dimerisation interface is formed by two identical surfaces, one on each subunit, coming together. This is how many proteins have evolved to be, and various functional advantages of homodimerisation (things like cooperativity, or higher stability, or easier folding, etc.) might have driven this evolution. But, generally speaking, for a feature to be selected for and optimised over time, it has to already exist in a primordial system, and perhaps if there is a bias towards that feature already at the start, it can be more easily picked up and optimised over time. Is there, therefore, an inherent bias for a random protein surface to associate with itself rather than with a different protein surface? This is a question posed in <a href="https://doi.org/10.1016/j.jmb.2006.11.020">this somewhat old paper</a> from Lukatsky et al. (the Shakhnovich lab). The paper contains a simplified model of a protein surface through which this question is investigated computationally (giving an affirmative answer), but there is also a common-sense argument, which I like, and which I think can be rephrased in the following way. What is special about a typical homodimerisation interface? That it has what crystallographers call  2-fold rotational symmetry (C2 symmetry), which means that all (or most) interactions are repeated twice. If residue X from subunit 1 interacts with residue Y from subunit 2, then residue Y from subunit 1 interacts with residue X from subunit 2. This repetition has implications for the likelihood of making homodimers, because, if you make certain favourable interactions by chance, they will be repeated twice. I think this is similar to what Lukatsky et al. say here: "it is more probable to symmetrically match a half of a random pattern with itself than with a different random pattern, which requires a full match". If we imagine that dimerisation is pattern-matching, it is more probable to match half of a pattern (& symmetry does the other half) than to match a full pattern (which is what has to happen during hetero-dimerisation)

I think that the repeptition idea also implies that the rate with which an interface becomes optimised over time is higher for a homodimer than a heterodimer - because a favourable mutation that strengthens the interface would have an amplified effect (would effectively be a double mutation). As Monod and his colleagues working on haemoglobin already observed, symmetry amplifies effects. Perhaps this is the common-sense rationale behind <a href="https://doi.org/10.1016/j.jmb.2009.10.044">another modelling study</a>, this time from G. E. Schulz and with a more sophisticated model tracking evolution over time. This second study concluded that a simulated evolutionary process that leads to homodimerisation has the probability that is by a factor of 100 or above higher than the probability of an analogous simulated process that leads to heterodimerisation.

