FROM ghcr.io/mspass-team/mspass:latest

LABEL maintainer="Ian Wang <yinzhi.wang.cug@gmail.com>"
ENV PFPATH=/test/pf

# Add cxx library
ADD cxx /parallel_pwmig/cxx
ADD data /parallel_pwmig/data
ADD setup.py /parallel_pwmig/setup.py
ADD python /parallel_pwmig/python
RUN MSPASS_HOME=/usr/local pip3 install /parallel_pwmig -v
RUN pip3 install vtk

