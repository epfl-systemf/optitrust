FROM ubuntu:22.04

RUN apt-get update && apt-get install -y opam clang clang-format libclang-dev llvm-dev libomp-dev pkg-config zlib1g-dev
# for C++ headers support:
RUN apt-get install -y libc++-dev


RUN opam init -y --disable-sandboxing
RUN opam switch -y create 4.14.1
RUN opam option depext-run-installs=true
RUN opam pin -y add menhirLib 20210419
RUN opam pin -y add pprint 20220103
RUN apt-get install -y autoconf libclang-cpp14-dev
RUN opam pin -y --assume-depexts add clangml 4.8.0
RUN opam install -y dune refl clangml pprint menhir menhirLib base64 ocamlbuild ocaml-lsp-server ppx_deriving
RUN ln -s /bin/true /bin/code