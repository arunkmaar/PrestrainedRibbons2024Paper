<TeXmacs|2.1.1>

<style|<tuple|tmarticle|preview-ref|smart-ref|comment|schola-font|number-long-article>>

<\body>
  <\hide-preamble>
    <assign|my-mult|*>

    <assign|t:|<syntax|:|\<cdot\>>>

    <assign|t.:|<syntax|<shift|\<therefore\>||.06em>|\<cdot\>>>

    \;
  </hide-preamble>

  <doc-data|<\doc-title>
    \ Solving for solutions in <with|font-family|tt|Auto-07p>
  </doc-title>|<doc-running-title|Wavelength selection in pre-strained
  ribbons>>

  <\doc-running-title>
    Solving for solutions in Auto-07p

    \;
  </doc-running-title>

  This file explains how to solve the boundary value problem (BVP) described
  in Section<nbsp>3.4 of the paper using <with|font-family|tt|Auto-07p><nbsp><cite|doedel2007auto>,
  a powerful differential equation solver that uses the arc length
  continuation.

  There are two cases discussed in Section 5 of the paper: the preferred long
  wavelength case and the short wavelength case. We have two folders
  corresponding to these cases in this directory:

  <\enumerate>
    <item><with|font-series|bold|Preferred Long Wavelength>

    <item><with|font-series|bold|Short Wavelength>
  </enumerate>

  In each folder, there are two sub-folders corresponding to respective
  solution scenarios:

  <\itemize>
    <item><with|font-series|bold|OnePerversionSolutions>

    <item><with|font-series|bold|TwelvePerversionSolutions>
  </itemize>

  These sub-folders contain the files required for solving the BVP. The
  computed solutions are shown in the videos in the parent
  <verbatim|Videos-ResultsFromNumericalSimulations> folder.

  <section*|Steps for obtaining solutions>

  The process of obtaining solutions involves two main steps:

  <\itemize>
    <item><with|font-series|bold|Solving the BVP in
    <with|font-family|tt|Auto-07p>>

    <item><with|font-series|bold|Post-processing the computed solutions>
  </itemize>

  We will explain these steps using the preferred long wavelength
  one-perversion case as an example. The procedure is similar for the other
  cases. Each subfolder contains <with|font-family|tt|Auto-07p> files and
  solutions are stored in the respective <with|font-family|tt|Output>.zip.
  You can use the long wavelength one-perversion case as a guide for
  post-processing solutions for other cases.

  <section|Solving the BVP in <with|font-family|tt|Auto-07p>>

  To solve the BVP, follow these instructions:

  <\itemize>
    <item>Go to the folder:

    \ <with|font-family|tt|/PreferredLongWavelength/OnePerversionSolutions>

    <item>Open <with|font-family|tt|Auto-07p> and run the
    <with|font-family|tt|solve.py> script \ by executing the command:

    <verbatim|execfile('solve.py')>
  </itemize>

  <\verbatim>
    \;

    \;
  </verbatim>

  The <with|font-family|tt|solve.py> script depends on two
  <with|font-family|tt|.f90> files, which set up the equilibrium equations,
  boundary conditions, and initial values for parameters. Specifically,
  <with|font-family|tt|Auto-07p> requires us to define the differential
  equations as a first-order system in the subroutine
  <with|font-family|tt|FUNC> in the <with|font-family|tt|.f90> files.

  The table below explains the variables used in the
  <with|font-family|tt|.f90> files and their corresponding quantities in the
  paper:

  <\equation*>
    <tabular*|<tformat|<cwith|2|7|1|1|cell-halign|c>|<cwith|2|7|3|3|cell-halign|c>|<cwith|1|1|1|1|cell-tborder|0ln>|<cwith|1|1|1|1|cell-bborder|1ln>|<cwith|2|2|1|1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|1|1|1|1|cell-rborder|0ln>|<cwith|1|1|2|2|cell-lborder|0ln>|<cwith|1|1|3|3|cell-tborder|0ln>|<cwith|1|1|3|3|cell-bborder|1ln>|<cwith|2|2|3|3|cell-tborder|1ln>|<cwith|1|1|3|3|cell-lborder|0ln>|<cwith|1|1|2|2|cell-rborder|0ln>|<cwith|1|1|3|3|cell-rborder|0ln>|<table|<row|<cell|<text|variable
    used in the <with|font-family|tt|.f90>
    files>>|<cell|>|<cell|<text|quantity in the
    paper>>>|<row|<cell|<text|<with|font-family|tt|theta>>>|<cell|>|<cell|\<theta\>>>|<row|<cell|<text|<with|font-family|tt|tau>><text|>>|<cell|>|<cell|\<tau\>>>|<row|<cell|<text|<with|font-family|tt|dTau>>>|<cell|>|<cell|<frac|\<mathd\>\<tau\>|\<mathd\>S>>>|<row|<cell|<text|<with|font-family|tt|d2Tau>>>|<cell|>|<cell|<frac|\<mathd\><rsup|2>\<tau\>|\<mathd\>S<rsup|2>>>>|<row|<cell|<text|<with|font-family|tt|d3TauByl>>>|<cell|>|<cell|<frac|1|\<ell\>><frac|\<mathd\><rsup|3>\<tau\>|\<mathd\>S<rsup|3>>>>|<row|<cell|<text|<with|font-family|tt|M>>>|<cell|>|<cell|M
    >>>>>
  </equation*>

  Also, we use variable <with|font-family|tt|Su<math|>><math|\<in\><around*|[|0,1|]>>
  as the independent variable for setting up the differential equations,
  where we define <with|font-family|tt|Su><math|=S/l>. The seventh
  differential equation in the subroutine <with|font-family|tt|FUNC> along
  with boundary seventh boundary condition in the subroutine
  <with|font-family|tt|BCND> specifies <with|font-family|tt|Su>.

  There are three parameters in the <with|font-family|tt|FUNC> subroutine:

  <\enumerate>
    <item><with|font-family|tt|epsilon> (<math|\<varepsilon\>> in the paper)
    controls the separation between ribbon ends.

    <item><with|font-family|tt|mExt> is the amplitude of the external
    sinusoidal moment used for the guess step explained below.

    <item>The third parameter, defined in subroutine
    <with|font-family|tt|PVLS>, measures <math|<around*|\<\|\|\>|\<tau\>|\<\|\|\>><rsub|<with|font-family|rm|rms>>>
    in the computed solutions. We refer to it as <with|font-family|tt|l2Tau>
    in the constants file <with|font-family|tt|c.guessStep> and the script
    file <with|font-family|tt|solve.py>.
  </enumerate>

  The first-order differential equations are expressed in the form
  <with|font-family|tt|F(i) = ...> in the <with|font-family|tt|.f90> files,
  where:

  <\equation*>
    <tabular*|<tformat|<cwith|1|-1|5|5|cell-halign|l>|<cwith|1|1|5|5|cell-tborder|0ln>|<cwith|1|1|5|5|cell-bborder|1ln>|<cwith|2|2|5|5|cell-tborder|1ln>|<cwith|1|1|5|5|cell-lborder|0ln>|<cwith|1|1|4|4|cell-rborder|0ln>|<cwith|1|1|5|5|cell-rborder|0ln>|<table|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|corresponding
    quantity in the paper>>|<row|<cell|<text|<with|font-family|tt|F(1)>>>|<cell|=>|<cell|<frac|\<mathd\>*<text|<with|font-family|tt|theta>>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|\<ell\>
    \<tau\>>>|<row|<cell|<text|<with|font-family|tt|F(2)>>>|<cell|=>|<cell|<frac|\<mathd\>*<with|font-family|tt|tau>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|\<ell\>*<frac|\<mathd\>*\<tau\>|\<mathd\>*S>>>|<row|<cell|<text|<with|font-family|tt|F(3)>>>|<cell|=>|<cell|<frac|\<mathd\>*<with|font-family|tt|dTau>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|<math|\<ell\>*<frac|\<mathd\><rsup|>*<rsup|2>\<tau\>|\<mathd\>*S<rsup|2>>>>>|<row|<cell|<text|<with|font-family|tt|F(4)>>>|<cell|=>|<cell|<frac|\<mathd\>*<text|<with|font-family|tt|d2Tau>>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|\<ell\>**<frac|\<mathd\><rsup|>*<rsup|3>\<tau\>|\<mathd\>*S<rsup|3>>>>|<row|<cell|<text|<with|font-family|tt|F(5)>>>|<cell|=>|<cell|<frac|\<mathd\>*<text|<with|font-family|tt|d3TauByl>>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|<frac|\<mathd\><rsup|>*<rsup|4>\<tau\>|\<mathd\>*S<rsup|4>>
    <around*|(|expression solved from eq. 3.10 in the
    paper|)>>>|<row|<cell|<text|<with|font-family|tt|F(6)>>>|<cell|=>|<cell|<frac|\<mathd\>*<text|<with|font-family|tt|M>>|\<mathd\>*<text|<with|font-family|tt|Su>>>>|<cell|=>|<cell|0
    <around*|(|the equilibrium equation|)>>>>>>
  </equation*>

  \;

  <with|font-series|bold|Loading procedure.> Using analytical calculations
  based on the dispersion relation in equation 4.5 of the paper, we
  anticipate that buckling points may be observed close to each other as we
  decrease the separation between ribbon ends. For example, in the preferred
  long wavelength case, the first twelve buckling points may occur between
  <math|\<varepsilon\>=-29.04> and <math|\<varepsilon\>=-31.02> as per Figure
  5.3(a).

  Owing to the difficulty in finding and selecting a solution branch
  corresponding to a specific number of perversions, we use the following
  procedure as echo in <with|font-family|tt|solving.py>:

  <\enumerate>
    <item><with|font-series|bold|Guess Step>: We first apply a sinusoidal
    external moment to create the desired number of perversions, in this case
    one perversion. This loading step, called <with|font-family|tt|guessStep>
    in <with|font-family|tt|solving.py>, increases the moment amplitude from
    0 to 100. Then, decrease the separation at a fixed moment until the
    ribbon ends is well below correpsonding buckling strain predicted for the
    desired mode. For example, decrease separation until
    <math|\<varepsilon\>=-45.9>. And then set the external moment to zero.

    <item><with|font-series|bold|Changing Separation>: With one perversion
    ribbon state at hand, we next perform a continuation in
    <math|\<varepsilon\>> to obtain the solution curve. We adjust constants
    and step sizes to ensure convergence. These values are specified in
    <with|font-family|tt|solve.py>. For explanations of various constants in
    <with|font-family|tt|solve.py> and <with|font-family|tt|c.guessStep>,
    refer to section A.2 of <cite|autoTutorial>.
  </enumerate>

  <section|Post-processing the computed solutions>

  After executing the script solve.py, we should see the following
  \ <with|font-family|tt|l2Tau> vs <with|font-family|tt|epsilon> plot on the
  screen.\ 

  <\big-figure|<image|l2TauvsEpsilonPlot.pdf|0.8par|||>>
    <with|font-family|tt|l2Tau> vs <with|font-family|tt|epsilon> plot that
    appears after the execution of script solve.py. \ 
  </big-figure>

  \;

  This plot corresponds to Figure 5.1 (with a negative sign for
  <math|\<varepsilon\>> in the paper). Towards the end of the end of the
  script solve.py, we save <with|font-family|tt|b.changeSeparation> file that
  stores the details of epsilon's increments used to compute the figure, and
  <with|font-family|tt|s.changeSeparation> that stores the solved ribbon
  states.

  The output from our <with|font-family|tt|Auto-07p> implementation is
  located inside the <with|font-family|tt|Output>.zip. You can unzip the
  folder and use the output if you want to circumvent solving the equations
  in <with|font-family|tt|Auto-07p> and want to directly use the solutions we
  got. For more details on how to interpret the columns in
  <with|font-family|tt|s.changeSeparation>, refer to section A.3 of
  <cite|autoTutorial>.

  For further post-processing, we use the Mathematica files in the
  <with|font-family|tt|Postprocessing> folder. We will use
  <with|font-family|tt|auto07pfilesParser.m> to import data from
  s.changeSeparation. We plot twist for a ribbon sate at
  <math|\<varepsilon\>=-29.07> in the Mathematica file.

  <\bibliography|bib|tm-plain|autoRib>
    <\bib-list|2>
      <bibitem*|1><label|bib-doedel2007auto>E.<nbsp>Doedel,
      A.<nbsp>R.<nbsp>Champneys, F.<nbsp>Dercole, T.<nbsp>F.<nbsp>Fairgrieve,
      Y.<nbsp>A.<nbsp>Kuznetsov, B.<nbsp>Oldeman, R.<nbsp>C.<nbsp>Paffenroth,
      B.<nbsp>Sandstede, X.<nbsp>J.<nbsp>Wang<localize|, and
      >C.<nbsp>H.<nbsp>Zhang. <newblock>Auto-07p: continuation and
      bifurcation software for ordinary differential equations.
      <newblock>2007.<newblock>

      <bibitem*|2><label|bib-autoTutorial>B.<nbsp>Sandstede<localize| and
      >D.<nbsp>Lloyd. <newblock>Using auto for stability problmes.
      <newblock><slink|https://github.com/sandstede-lab/Auto07p/blob/master/auto07p_tutorial_spatial_pattern_formation/auto07p_tutorial_spatial_pattern_formation.pdf>.
      <newblock>Accessed: 2024-08-15.<newblock>
    </bib-list>
  </bibliography>

  \;
</body>

<\initial>
  <\collection>
    <associate|font-base-size|8>
    <associate|info-flag|detailed>
    <associate|math-font|math-schola>
    <associate|page-medium|paper>
    <associate|preamble|false>
  </collection>
</initial>

<\attachments>
  <\collection>
    <\associate|bib-bibliography>
      <\db-entry|+Nj1ldIT1klNWWhg|article|huang2018differential>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1704363619>
      <|db-entry>
        <db-field|author|Changjin <name|Huang><name-sep>Zilu
        <name|Wang><name-sep>David <name|Quinn><name-sep>Subra
        <name|Suresh><name-sep>K Jimmy <name|Hsia>>

        <db-field|title|Differential growth and shape formation in plant
        organs>

        <db-field|journal|Proceedings of the National Academy of Sciences>

        <db-field|year|2018>

        <db-field|volume|115>

        <db-field|number|49>

        <db-field|pages|12359\U12364>
      </db-entry>

      <\db-entry|+1g0I8RJKyWdEmmS|article|savin2011growth>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1712823672>
      <|db-entry>
        <db-field|author|Thierry <name|Savin><name-sep>Natasza A
        <name|Kurpios><name-sep>Amy E <name|Shyer><name-sep>Patricia
        <name|Florescu><name-sep>Haiyi <name|Liang><name-sep>L
        <name|Mahadevan><name-sep>Clifford J <name|Tabin>>

        <db-field|title|On the growth and form of the gut>

        <db-field|journal|Nature>

        <db-field|year|2011>

        <db-field|volume|476>

        <db-field|number|7358>

        <db-field|pages|57\U62>

        <db-field|doi|10.1038/nature10277>

        <db-field|publisher|Nature Publishing Group UK London>
      </db-entry>

      <\db-entry|+1g0I8RJKyWdEmmT|misc|livescienceOctopusesChange>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1712823672>
      <|db-entry>
        <db-field|author|Harry <name|Baker>>

        <db-field|title|How do octopuses change color?>

        <db-field|howpublished|<slink|https://www.livescience.com/how-do-octopuses-change-color>>

        <db-field|year|2022>

        <db-field|note|[Accessed 05-04-2024]>

        <db-field|url|<slink|https://www.livescience.com/how-do-octopuses-change-color>>
      </db-entry>

      <\db-entry|+21Pg386qSekctPO|article|huang2012spontaneous>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1711710916>
      <|db-entry>
        <db-field|author|Jiangshui <name|Huang><name-sep>Jia
        <name|Liu><name-sep>Benedikt <name|Kroll><name-sep>Katia
        <name|Bertoldi><name-sep>David R <name|Clarke>>

        <db-field|title|Spontaneous and deterministic three-dimensional
        curling of pre-strained elastomeric bi-strips>

        <db-field|journal|Soft Matter>

        <db-field|year|2012>

        <db-field|volume|8>

        <db-field|number|23>

        <db-field|pages|6291\U6300>

        <db-field|doi|10.1039/C2SM25278C>

        <db-field|publisher|Royal Society of Chemistry>

        <db-field|url|<slink|http://dx.doi.org/10.1039/C2SM25278C>>
      </db-entry>

      <\db-entry|+21Pg386qSekctPN|article|liu2014structural>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1711710803>
      <|db-entry>
        <db-field|author|Jia <name|Liu><name-sep>Jiangshui
        <name|Huang><name-sep>Tianxiang <name|Su><name-sep>Katia
        <name|Bertoldi><name-sep>David R <name|Clarke>>

        <db-field|title|Structural transition from helices to hemihelices>

        <db-field|journal|PloS one>

        <db-field|year|2014>

        <db-field|volume|9>

        <db-field|number|4>

        <db-field|pages|1\U7>

        <db-field|month|04>

        <db-field|doi|10.1371/journal.pone.0093183>

        <db-field|publisher|Public Library of Science San Francisco, USA>

        <db-field|url|<slink|https://doi.org/10.1371/journal.pone.0093183>>
      </db-entry>

      <\db-entry|+17jEphbnkOHqlwr|article|gomez2023twisting>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1705575818>
      <|db-entry>
        <db-field|author|Michael <name|Gomez><name-sep>Pedro M.
        <name|Reis><name-sep>Basile <name|Audoly>>

        <db-field|title|Twisting instabilities in elastic ribbons with
        inhomogeneous pre-stress: a macroscopic analog of thermodynamic phase
        transition>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2023>

        <db-field|volume|181>

        <db-field|pages|105420>

        <db-field|doi|<slink|https://doi.org/10.1016/j.jmps.2023.105420>>

        <db-field|issn|0022-5096>
      </db-entry>

      <\db-entry|+21Pg386qSekctPQ|article|lestringant2017elastic>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1711715036>
      <|db-entry>
        <db-field|author|Claire <name|Lestringant><name-sep>Basile
        <name|Audoly>>

        <db-field|title|Elastic rods with incompatible strain: macroscopic
        versus microscopic buckling>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2017>

        <db-field|volume|103>

        <db-field|pages|40\U71>

        <db-field|publisher|Elsevier>
      </db-entry>

      <\db-entry|+Nj1ldIT1klNWWhi|article|audoly2021one>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1704363619>
      <|db-entry>
        <db-field|author|B. <name|Audoly><name-sep>S. <name|Neukirch>>

        <db-field|title|A one-dimensional model for elastic ribbons: a little
        stretching makes a big difference>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2021>

        <db-field|volume|153>

        <db-field|pages|104457>
      </db-entry>

      <\db-entry|+Nj1ldIT1klNWWhf|article|lestringant2020asymptotically>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1704363619>
      <|db-entry>
        <db-field|author|C. <name|Lestringant><name-sep>B. <name|Audoly>>

        <db-field|title|Asymptotically exact strain-gradient models for
        nonlinear slender elastic structures: a systematic derivation method>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2020>

        <db-field|volume|136>

        <db-field|pages|103730>
      </db-entry>

      <\db-entry|+21Pg386qSekctPP|article|bardenhagen1994derivation>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1711715036>
      <|db-entry>
        <db-field|author|Scott <name|Bardenhagen><name-sep>Nicolas
        <name|Triantafyllidis>>

        <db-field|title|Derivation of higher order gradient continuum
        theories in 2, 3-d non-linear elasticity from periodic lattice
        models>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|1994>

        <db-field|volume|42>

        <db-field|number|1>

        <db-field|pages|111\U139>

        <db-field|doi|<slink|https://doi.org/10.1016/0022-5096(94)90051-5>>

        <db-field|publisher|Elsevier>

        <db-field|url|<slink|https://www.sciencedirect.com/science/article/pii/0022509694900515>>
      </db-entry>

      <\db-entry|+2URFiYCV2pgZIEg|article|durand2022predictive>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1712249062>
      <|db-entry>
        <db-field|author|Baptiste <name|Durand><name-sep>Arthur
        <name|Lebée><name-sep>Pierre <name|Seppecher><name-sep>Karam
        <name|Sab>>

        <db-field|title|Predictive strain-gradient homogenization of a
        pantographic material with compliant junctions>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2022>

        <db-field|volume|160>

        <db-field|pages|104773>

        <db-field|doi|<slink|https://doi.org/10.1016/j.jmps.2021.104773>>

        <db-field|publisher|Elsevier>
      </db-entry>

      <\db-entry|+Nj1ldIT1klNWWhe|article|kumar2023asymptotic>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1704363619>
      <|db-entry>
        <db-field|author|A. <name|Kumar><name-sep>B.
        <name|Audoly><name-sep>C. <name|Lestringant>>

        <db-field|title|Asymptotic derivation of a higher-order
        one-dimensional model for tape springs>

        <db-field|journal|Philosophical Transactions of the Royal Society A>

        <db-field|year|2023>

        <db-field|volume|381>

        <db-field|number|2244>

        <db-field|pages|20220028>
      </db-entry>

      <\db-entry|+1g0I8RJKyWdEmmR|article|koiter1960consistent>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1712823672>
      <|db-entry>
        <db-field|author|WT <name|Koiter>>

        <db-field|title|A consistent first approximation in the general
        theory of thin elastic shells>

        <db-field|journal|The theory of thin elastic shells>

        <db-field|year|1960>

        <db-field|pages|12\U33>

        <db-field|publisher|North-Holland Amsterdam>
      </db-entry>

      <\db-entry|+V99NGjo2EXdXlAT|article|audoly2016analysis>
        <db-field|contributor|arunkumar>

        <db-field|modus|imported>

        <db-field|date|1704938181>
      <|db-entry>
        <db-field|author|B. <name|Audoly><name-sep>J. W. <name|Hutchinson>>

        <db-field|title|Analysis of necking based on a one-dimensional model>

        <db-field|journal|Journal of the Mechanics and Physics of Solids>

        <db-field|year|2016>

        <db-field|volume|97>

        <db-field|pages|68\U91>

        <db-field|doi|<slink|https://doi.org/10.1016/j.jmps.2015.12.018>>

        <db-field|publisher|Elsevier>
      </db-entry>
    </associate>
  </collection>
</attachments>

<\references>
  <\collection>
    <associate|auto-1|<tuple|<with|mode|<quote|math>|\<bullet\>>|1>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|2.1|2>>
    <associate|auto-5|<tuple|2.1|3>>
    <associate|bib-autoTutorial|<tuple|2|3>>
    <associate|bib-doedel2007auto|<tuple|1|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      doedel2007auto

      autoTutorial

      autoTutorial
    </associate>
    <\associate|figure>
      <tuple|normal|<\surround|<hidden-binding|<tuple>|2.1>|>
        <with|font-family|<quote|tt>|l2Tau> vs
        <with|font-family|<quote|tt>|epsilon> plot that appears after the
        execution of script solve.py. \ 
      </surround>|<pageref|auto-4>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-shape|<quote|small-caps>|Steps
      for obtaining solutions> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-shape|<quote|small-caps>|1.<space|2spc>Solving
      the BVP in <with|font-family|<quote|tt>|Auto-07p>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-shape|<quote|small-caps>|2.<space|2spc>Post-processing
      the computed solutions> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-shape|<quote|small-caps>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <pageref|auto-5><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>