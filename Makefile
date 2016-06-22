#xfig -specialtext -latexfonts -startlatexFont default
#export "Combined PDF/LaTeX (both parts)"
%.pdf %.pdf_t: %.fig
	fig2dev -L pdftex $*.fig $*.pdf
	fig2dev -L pdftex_t -p $*.pdf $*.fig $*.pdf_t
%.png: %.gpt
	gnuplot $*.gpt
%.pdf_t: %.eps %.tex # Transform from gnuplot with terminal=epslatex to pdf
	epstopdf $*.eps
	cp $*.tex $*.pdf_t
%.pdf: %.tex
	rubber --pdf $<
AUGFIGS= rho-tilde.png u-tilde.png rate.png rho-hat.png rate-hat.png count.png back-proj.png augfig2.pdf augfig1.pdf

aug07.pdf: aug07.tex $(AUGFIGS)

TMT.pdf: TMT.tex TMTF.pdf_t TMTF.pdf V_0.pdf U_0.pdf R_true.pdf Count.pdf R_hat.pdf U_hat.pdf V_hat.pdf

ToDo9_07.pdf: ToDo9_07.tex rho-tilde.png rho-hat.png back-proj.png overview.pdf_t 
Plan_08.pdf: Plan_08.tex rho-tilde.png rho-hat.png back-proj.png overview.pdf_t 

fuse_scheme: # See ~/fuse/README
	mkdir -p fuse_scheme
fuse_scheme/A_LUISbrick_runlist.txt: fuse_scheme
	sshfs -o idmap=user afraser@yks.lanl.gov:/n/toaster/u/agreen/luis_overhead_data fuse_scheme
amf.so: amf.c
	python setup.py build
# Rule for simulations
%true.eps %true.tex %hat.eps %hat.tex: %S.fig_com GUILM.py
	python GUILM.py simulation $<
# Rule for experimental data
%hat.eps %hat.tex: %E.fig_com GUILM.py
	python GUILM.py data $<
LMTMT.pdf: LMTMT.tex exp1hat.pdf_t sim1hat.pdf_t sim1true.pdf_t
	pdflatex $<

RELEASE_FILES = Makefile geometryLM.py GUILM.py setup.py amf.c	\
LMTMT.tex README LMTMT.pdf guide.tex guide.pdf #tmp/AG_rotated.mat

lmt4luis.tar.bz2: ${RELEASE_FILES}
	tar -cjf $@ $^

DERIVED = $(patsubst %, *%, pdf log aux bbl blg snm nav toc pdf_t eps pyc ~\
.so png .out)
.PHONY: clean
clean:
	rm ${DERIVED}
	rm -rf build fuse_scheme

# Local Variables:
# mode: makefile
# End:
