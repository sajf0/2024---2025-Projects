\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ArabicClass}[2024/02/09 Example LaTeX class]

\LoadClass[14pt]{extreport}

\RequirePackage{tasks,enumitem}
\RequirePackage[english]{babel}
\RequirePackage[a4paper, width=150mm,top=35mm,bottom=35mm]{geometry}
\RequirePackage{bm,fixmath} 
\RequirePackage{amsthm,amsfonts,amsmath,amssymb,mathtools,nicefrac,xfrac,siunitx,nccmath}
\RequirePackage{titletoc,titlesec,fmtcount}

\RequirePackage[nodisplayskipstretch]{setspace}
\RequirePackage{graphicx,xcolor,hyperref,float}
\RequirePackage{tcolorbox}
\RequirePackage{pgfplots}
\pgfplotsset{compat=1.18}
\RequirePackage{titling,authblk,chngcntr,subcaption}
\RequirePackage[twoside]{fancyhdr}
\setlength{\headheight}{16.7pt}
\pagestyle{fancy}
\fancyhead{}
\fancyfoot{}
\fancyhead[LO,RE]{\thepage}
\fancyhead[CE]{\small \leftmark}
\fancyhead[CO]{\small \thechapter}

\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{#1}}{}}

\RequirePackage{polyglossia}
\setmainlanguage[numerals=maghrib]{arabic}
\setotherlanguage{english}
\newfontfamily\arabicfont[Script=Arabic]{Times New Roman}
\usepackage[T1]{fontenc}
\usepackage{times}
\usepackage[lite, subscriptcorrection, zswash]{mtpro2}

\newcommand{\en}{\textenglish}
\newcommand{\ar}{\textarabic}
\newcommand{\maghrib}{\setmainlanguage[numerals=maghrib]{arabic}}
\newcommand{\mashriq}{\setmainlanguage[numerals=mashriq]{arabic}}
\rightfootnoterule

\RequirePackage{biblatex}

\titleformat{\chapter}[display]
{\normalfont\Huge\bfseries\centering}
  {\thechapter}
  {10pt}
{\Huge\bfseries\centering}

\titlespacing*{\chapter}{0pt}{150pt}{150pt}
\titleclass{\subchapter}{top}[\chapter]
\newcounter{subchapter}
\titleformat{\subchapter}[display]
{\normalfont\Large\bfseries\centering}
  {\thechapter}
  {10pt}
{\Large\bfseries\centering}
\titlespacing*{\subchapter}{0pt}{0pt}{20pt}

\titleformat{\section}[block]{\color{red}\normalfont\bfseries\raggedleft}{\thesection}{8pt}{}
%\titlespacing*{\section}{5pt}{\parskip}{0pt}
\titleformat{\subsection}[block]{\normalfont\bfseries\raggedleft}{\thesubsection}{5pt}{}
%\titlespacing*{\subsection}{5pt}{\parskip}{0pt}

\renewcommand{\thetable}{\arabic{chapter}-\arabic{table}}
\renewcommand{\thefigure}{\arabic{chapter}-\arabic{figure}}

\renewcommand{\thechapter}  {\ifnum\value{chapter}=1 الفصل الأول
   \else
     \ifnum\value{chapter}=2 الفصل الثاني
      \else
      \ifnum\value{chapter}=3 الفصل الثالث
      \else
      \ifnum\value{chapter}=4 الفصل الرابع
      \else
      \ifnum\value{chapter}=5 الفصل الخامس
      \else
      \ifnum\value{chapter}=6 الفصل السادس
      \else
      \ifnum\value{chapter}=7 الفصل السابع
      \else
     \Numberstring{chapter}
      \fi
      \fi
      \fi
      \fi
      \fi
      \fi
   \fi}
   
\renewcommand\thesection{\arabic{chapter}-\arabic{section}}
\renewcommand\thesubsection{\arabic{chapter}-\arabic{section}-\arabic{subsection}}

\renewcommand{\tableofcontents}{%
    \begingroup
    \let\clearpage\relax 
    \subchapter*{المحتويات}
    \@starttoc{toc}
    \endgroup
}

\titlecontents{chapter}[0em]{\bfseries\vspace{15pt}}{\thecontentslabel: }{}{\hfill}

\titlecontents{chapter*}[0em]{\bfseries\vspace{15pt}}{}{}{\hfill\normalfont\contentspage}

\titlecontents{section}[1em]
{\normalfont}
{\thecontentslabel\hspace*{5pt} }
{}{\titlerule*[0.5pc]{.}\contentspage}

\hypersetup{
colorlinks=true,
linkcolor=black,
citecolor=red,
}

\renewcommand*\newunitpunct{\addcomma\space}
\makeatletter
\def\abx@comma{\ifdim\lastkern>\z@\unkern\fi\abx@puncthook{،}}

\DeclareMathOperator{\K}{\,K^o}
\DeclareMathOperator*{\Res}{Res\,}
\DeclareMathOperator{\E}{E\,}
\DeclareMathOperator{\Var}{Var\,}
\DeclareMathOperator{\pr}{P_r\,}

\newcommand{\der}[2]{\frac{d#1}{d#2}}
\newcommand{\pder}[2]{\dfrac{\partial#1}{\partial#2}}

\newcommand{\LL}{\mathcal{L}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\boldN}{\pmb{\mathbb{N}}}
\newcommand{\boldZ}{\pmb{\mathbb{Z}}}
\newcommand{\boldQ}{\pmb{\mathbb{Q}}}
\newcommand{\boldR}{\pmb{\mathbb{R}}}
\newcommand{\boldC}{\pmb{\mathbb{C}}}

\newcommand{\pbracket}[1]{\left(#1\right)}
\newcommand{\sbracket}[1]{\left[#1\right]}
\newcommand{\cbracket}[1]{\left\{#1\right\}}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\subs}[1]{\left.#1\right|}

\newtheoremstyle{mystyle}
  {\topsep}
  {\topsep}
  {\color{black}}
  {}
  {\bfseries}
  {\vspace{5pt}\\}
  {.5ex}  
  {#1 #2 \textbf{#3}}
  
\theoremstyle{mystyle}
\newtheorem{definition}{تعريف}[section]

\newtheorem{note}[definition]{ملاحظة}
\newtheorem{example}[definition]{مثال}
\newtheorem*{solution}{الحل}
\newtheorem{theorem}[definition]{مبرهنة}
\newenvironment{myproof}{\noindent\textbf{البرهان}\vspace{5pt}\\}{ $\qedsymbol$}

\renewcommand{\thenote}{\arabic{note}}
\renewcommand{\thedefinition}{\arabic{chapter}-\arabic{section}-\arabic{definition}}
%\renewcommand{\theexample}{\arabic{chapter}-\arabic{example}}
%\renewcommand{\thetheorem}{\arabic{chapter}-\arabic{theorem}}


\endinput