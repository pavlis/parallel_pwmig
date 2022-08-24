FROM ghcr.io/mspass-team/mspass:latest

LABEL maintainer="Ian Wang <yinzhi.wang.cug@gmail.com>"

# Add cxx library
ADD cxx /parallel_pwmig/cxx
RUN cd /parallel_pwmig/cxx \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && make install \ 
    && rm -rf ../build

ADD setup.py /parallel_pwmig/setup.py
ADD python /parallel_pwmig/python
RUN pip3 install /parallel_pwmig -v

