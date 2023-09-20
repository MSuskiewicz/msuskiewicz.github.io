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
A lot of proteins we see in nature make homodimers in which the dimerisation interface is formed by two identical surfaces, one on each subunit, coming together. This is how many proteins have evolved to be, and various functional advantages of homodimerisation might have driven this evolution. But, generally speaking, for a feature to be selected for and optimised over time, it has to already exist in a primordial system, and perhaps if there is a bias towards that feature already at the start, it can be more easily picked up and optimised over time. Is there, therefore, an inherent bias for a random protein surface to associate with itself rather than with a different protein surface? This is a question posed in <a href="https://doi.org/10.1016/j.jmb.2006.11.020">this somewhat old paper</a> from Lukatsky et al. (the Shakhnovich lab). The paper contains a simplified model of a protein surface through which this question is investigated (giving an affirmative answer), but there is also a common-sense argument, which I like, and which I think can be rephrased in the following way. What is special about a typical homodimerisation interface? That all interactions are repeated twice. If residue X from subunit 1 interacts with residue Y from subunit 2, then residue Y from subunit 1 interacts with residue X from subunit 2. This repetition has implications for the likelihood of making homodimers, because, if you by chance make certain favourable interactions, they will be repeated twice. I think this is similar to what Lukatsky et al. say here: "it is more probable to symmetrically match a half of a random pattern with itself than with a different random pattern, which requires a full match". I think that the repeptition idea also implies that an effect of a single point mutations at a homodimerisation interface is stronger than at a heterodimerisation interface, because in a homodimer a single mutation is effectively a double mutation. If a single mutation strengthens the dimerisation, this effect will be stronger for a homodimer than a heterodimer - hence not only are wvery weak homodimers more likely by chance than very weak heterodimers, but also homodimers might be strengthened faster during evolution than heterodimers are. Perhaps this is the common-sense rationale behind <a href="https://doi.org/10.1016/j.jmb.2009.10.044">another modelling study</a>, this time with a more sophisticated model tracking evolution over time, that concluded that a simulated evolutionary process that leads to homodimerisation has the probability that is by a factor of 100 or above higher than the probability of a process that leads to heterodimerisation.

