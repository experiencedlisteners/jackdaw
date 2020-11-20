FROM rigetti/lisp
RUN apt-get update && apt-get -y install emacs vim
RUN echo '(ql:quickload "quicklisp-slime-helper")' | sbcl

COPY emacs /root/.emacs
# Replace with pulls from repos
COPY ./jackdaw /jackdaw 
COPY ./idyom /idyom 

WORKDIR /jackdaw/jackdaw/
ENV JACKDAW_ROOT=/jackdaw/
ENV IDYOM_ROOT=/idyom/
# Trigger ql dependency installation
RUN /jackdaw/jackdaw/cli.lisp
VOLUME ./data /data
ENTRYPOINT ["/jackdaw/jackdaw/cli.lisp"]
#ENTRYPOINT ["cli.lisp"]
# Install slime through quicklisp
# Add quicklisp lines to .emacs
# Mount models directory

