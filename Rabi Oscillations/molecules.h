#ifndef MOLECULES_H
#define MOLECULES_H


struct mol {
        int molname;
        double pos[3];
        double vel[3];
        double tinit;
        int state;
};

#endif
