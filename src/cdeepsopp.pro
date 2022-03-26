TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += "../wup/cpp/include"

SOURCES += \
        main.cpp

HEADERS += \
    cdeepso.hpp \
    cdeepso_params.hpp \
    functions.hpp \
    ntuplecdeepso.hpp \
    operations.hpp \
    population.hpp \
    utils.hpp \
    weight.hpp
