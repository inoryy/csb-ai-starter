#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("omit-frame-pointer")
#pragma GCC optimize("unroll-loops")

#include "stdio.h"
#include "math.h"
#include <iostream>
#include <algorithm>
#include <memory>
#include <chrono>
#include <vector>

using namespace std;
using namespace std::chrono;

high_resolution_clock::time_point now = high_resolution_clock::now();
#define TIME duration_cast<duration<double>>(high_resolution_clock::now() - now).count()

class Point;
class Unit;
class Pod;
class Collision;
class Checkpoint;
class Solution;
class Bot;
void load();
void play();
void print_move(int, float, Pod*);

constexpr int CP  = 0;
constexpr int POD = 1;
constexpr int DEPTH = 6;
constexpr float SHIELD_PROB = 10;
constexpr int MAX_THRUST = 100;

constexpr float E = 0.00001;

int r = -1;
int turn = 0;
int sols_ct = 0;
bool is_p2 = false;
int cp_ct, laps;

Pod* pods[4];
Checkpoint* cps[10];

inline int fastrand() {
    static unsigned int g_seed = 42;
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

inline int rnd(int b) {
    return fastrand() % b;
}

inline int rnd(int a, int b) {
    return a + rnd(b - a + 1);
}

class Collision {
public:
    Unit* a;
    Unit* b;
    float t;

    Collision() {}

    Collision(Unit* a, Unit* b, float t) {
        this->a = a;
        this->b = b;
        this->t = t;
    }
};

class Point {
public:
    float x, y;

    Point() {};

    Point(float x, float y) {
        this->x = x;
        this->y = y;
    }

    inline virtual float dist(Point p) {
        return sqrt(pow((x - p.x), 2) + pow((y - p.y), 2));
    }

    inline virtual float dist(Point* p) {
        return sqrt(dist2(p));
    }

    inline virtual float dist2(Point* p) {
        return pow((x - p->x), 2) + pow((y - p->y), 2);
    }

    Point closest(Point* a, Point* b) {
        float da = b->y - a->y;
        float db = a->x - b->x;
        float c1 = da*a->x + db*a->y;
        float c2 = -db*x + da*y;
        float det = da*da + db*db;

        float cx, cy;
        if (det != 0) {
            cx = (da*c1 - db*c2) / det;
            cy = (da*c2 + db*c1) / det;
        } else {
            cx = x, cy = y;
        }

        return Point(cx, cy);
    }
};

class Unit: public Point {
private:
    float cache[5];

public:
    int id, type;
    float r, vx, vy;

    virtual void bounce(Unit* u) {};

    inline float collision_time(Unit* u) {
        if (vx == u->vx && vy == u->vy) {
            return -1;
        }

        float sr2 = u->type == CP ? 357604 : 640000;

        float dx = x - u->x;
        float dy = y - u->y;
        float dvx = vx - u->vx;
        float dvy = vy - u->vy;
        float a = dvx*dvx + dvy*dvy;

        if (a < E) return -1;

        float b = -2.0*(dx*dvx + dy*dvy);
        float delta = b*b - 4.0*a*(dx*dx + dy*dy - sr2);

        if (delta < 0.0) return -1;

        float t = (b - sqrt(delta))*(1.0/(2.0*a));

        if (t <= 0.0 || t > 1.0) return -1;

        return t;
    }

    void save() {
        cache[0] = x;
        cache[1] = y;
        cache[2] = vx;
        cache[3] = vy;
    }

    void load() {
        x = cache[0];
        y = cache[1];
        vx = cache[2];
        vy = cache[3];
    }
};

class Checkpoint: public Unit {
public:
    Checkpoint(int id, float x, float y) {
        this->id = id;
        this->x = x;
        this->y = y;

        this->vx = this->vy = 0;
        this->type = CP;
        this->r = 600;
    }

    void bounce(Unit*) {}
};

class Pod: public Unit {
public:
    float angle = -1;
    float next_angle = -1;
    bool has_boost;
    int ncpid, checked, timeout, shield;
    Pod* partner;

    // TODO maybe replace cache array with primitives?
    float cache[10];

    Pod(int id) {
        this->id = id;
        this->r = 400;
        this->type = POD;
        this->ncpid = 1;
        // TODO move timeout to global/team var
        this->timeout = 100;
        this->has_boost = true;
        this->checked = this->shield = 0;
    }

    float score() {
        return checked*50000 - this->dist(cps[this->ncpid]);
    }

    void apply(int thrust, float angle) {
        angle = max((float)-18., min((float)18., angle));
        this->angle += angle;
        if (this->angle >= 360.) {
            this->angle = this->angle - 360.;
        } else if (this->angle < 0.0) {
            this->angle += 360.;
        }

        if (thrust == -1) {
            this->shield = 4;
        } else {
            boost(thrust);
        }
    }

    void rotate(Point* p) {
        float a = diff_angle(p);
        a = max((float)-18., min((float)18., a));

        angle += a;
        if (angle >= 360.) {
            angle = angle - 360.;
        } else if (angle < 0.0) {
            angle += 360.;
        }
    }

    void boost(int thrust) {
        if (shield > 0) return;

        float ra = angle * M_PI / 180.0;

        vx += cos(ra) * thrust;
        vy += sin(ra) * thrust;
    }

    void move(float t) {
        x += vx * t;
        y += vy * t;
    }

    void end() {
        x = round(x);
        y = round(y);
        vx = trunc(vx * 0.85);
        vy = trunc(vy * 0.85);

        if (checked >= cp_ct * laps) {
            ncpid = 0;
            checked = cp_ct * laps;
        }
        timeout--;
        if (shield > 0) shield--;
    }

    void bounce(Unit* u) {
        if (u->type == CP) {
            checked += 1;
            timeout = partner->timeout = 100;
            ncpid = (ncpid + 1) % cp_ct;
            return;
        }

        bounce_w_pod(static_cast<Pod*>(u));
    }

    void bounce_w_pod(Pod* u) {
        float m1 = shield == 4 ? 10. : 1.;
        float m2 = u->shield == 4 ? 10. : 1.;
        float mcoeff = (m1 + m2) / (m1 * m2);

        float nx = x - u->x;
        float ny = y - u->y;
        float dst2 = nx*nx + ny*ny;
        float dvx = vx - u->vx;
        float dvy = vy - u->vy;
        float prod = (nx*dvx + ny*dvy) / (dst2 * mcoeff);
        float fx = nx * prod;
        float fy = ny * prod;
        float m1_inv = 1.0 / m1;
        float m2_inv = 1.0 / m2;

        vx -= fx * m1_inv;
        vy -= fy * m1_inv;
        u->vx += fx * m2_inv;
        u->vy += fy * m2_inv;

        float impulse = sqrt(fx*fx + fy*fy);
        if (impulse < 120.) {
            float df = 120.0 / impulse;
            fx *= df;
            fy *= df;
        }

        vx -= fx * m1_inv;
        vy -= fy * m1_inv;
        u->vx += fx * m2_inv;
        u->vy += fy * m2_inv;
    }

    inline float diff_angle(Point* p) {
        float a = get_angle(p);
        float right = angle <= a ? a - angle : 360. - angle + a;
        float left = angle >= a ? angle - a : angle + 360. - a;

        if (right < left) {
            return right;
        }

        return -left;
    }

    inline float get_angle(Point* p) {
        float d = this->dist(p);
        float dx = (p->x - x) / d;
        float dy = (p->y - y) / d;

        float a = acos(dx) * 180 / M_PI;

        if (dy < 0) {
            a = 360 - a;
        }

        return a;
    }

    void update(int x, int y, int vx, int vy, float angle, int ncpid) {
        if (shield > 0) shield--;
        if (ncpid != this->ncpid) {
            timeout = partner->timeout = 100;
            checked++;
        } else {
            timeout--;
        }

        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
        this->ncpid = ncpid;

        if (is_p2 && id > 1) swap(angle, this->next_angle);
        this->angle = angle;
        if (::r == 0) this->angle = 1 + diff_angle(cps[1]);
        save();
    }

    void update(int shield, bool has_boost) {
        this->shield = shield;
        this->has_boost = has_boost;
    }

    void save() {
        Unit::save();
        cache[0] = ncpid;
        cache[1] = checked;
        cache[2] = timeout;
        cache[3] = shield;
        cache[4] = angle;
        cache[5] = has_boost;
    }

    void load() {
        Unit::load();
        ncpid   = cache[0];
        checked = cache[1];
        timeout = cache[2];
        shield  = cache[3];
        angle   = cache[4];
        has_boost = cache[5];
    }
};

class Solution {
public:
    float score = -1;
    int thrusts[DEPTH*2];
    float angles[DEPTH*2];

    Solution(bool with_rnd = false) {
        if (with_rnd) randomize();
    }

    void shift() {
        for (int i = 1; i < DEPTH; i++) {
            angles[i-1]        = angles[i];
            thrusts[i-1]       = thrusts[i];
            angles[i-1+DEPTH]  = angles[i+DEPTH];
            thrusts[i-1+DEPTH] = thrusts[i+DEPTH];
        }
        randomize(DEPTH-1, true);
        randomize(2*DEPTH-1, true);
        score = -1;
    }

    void mutate() {
        randomize(rnd(2*DEPTH));
    }

    void mutate(Solution* child) {
        copy(begin(angles), end(angles), begin(child->angles));
        copy(begin(thrusts), end(thrusts), begin(child->thrusts));

        child->mutate();
        child->score = -1;
    }

    void randomize(int idx, bool full = false) {
        int r = rnd(2);
        if (full || r == 0) angles[idx] = max(-18, min(18, rnd(-40, 40)));

        if (full || r == 1) {
            if (rnd(100) >= SHIELD_PROB) {
                thrusts[idx] = max(0, min(MAX_THRUST, rnd((int) -0.5*MAX_THRUST, 2*MAX_THRUST)));
            } else {
                thrusts[idx] = -1;
            }
        }
        score = -1;
    }

    void randomize() {
        for (int i = 0; i < 2*DEPTH; i++) randomize(i, true);
    }
};

class Bot {
public:
    int id = 0;

    Bot() {};

    Bot(int id) {
        this->id = id;
    }

    virtual void move() = 0;

    Pod* runner() {
        return runner(pods[id], pods[id+1]);
    }

    Pod* blocker() {
        return blocker(pods[id], pods[id+1]);
    }

    Pod* runner(Pod* pod0, Pod* pod1) {
        return pod0->score() - pod1->score() >= -1000 ? pod0 : pod1;
    }

    Pod* blocker(Pod* pod0, Pod* pod1) {
        return runner(pod0, pod1)->partner;
    }
};

class ReflexBot : public Bot {
public:
    ReflexBot() {}

    ReflexBot(int id) {
        this->id = id;
    }

    void move() {
        move_runner();
        move_blocker();
    }

    void move_as_main() {
        move_runner(true);
        move_blocker(true);
    }

    void move_runner(bool for_output = false) {
        Pod* pod = !for_output ? runner() : pods[0];

        Checkpoint* cp = cps[pod->ncpid];
        Point t(cp->x - 3*pod->vx, cp->y - 3*pod->vy);
        float raw_angle = pod->diff_angle(&t);

        int thrust = abs(raw_angle) < 90 ? MAX_THRUST : 0;
        float angle = max((float) -18, min((float) 18, raw_angle));

        if (!for_output) pod->apply(thrust, angle);
        else print_move(thrust, angle, pod);
    }

    void move_blocker(bool for_output = false) {
        Pod* pod = !for_output ? blocker() : pods[1];

        Checkpoint* cp = cps[pod->ncpid];
        Point t(cp->x - 3*pod->vx, cp->y - 3*pod->vy);
        float raw_angle = pod->diff_angle(&t);

        int thrust = abs(raw_angle) < 90 ? MAX_THRUST : 0;
        float angle = max((float) -18, min((float) 18, raw_angle));

        if (!for_output) pod->apply(thrust, angle);
        else print_move(thrust, angle, pod);
    }
};

class SearchBot : public Bot {
public:
    Solution sol;
    vector<Bot*> oppBots;

    SearchBot() {}

    SearchBot(int id) {
        this->id = id;
    }

    void move(Solution* sol) {
        pods[id]->apply(sol->thrusts[turn], sol->angles[turn]);
        pods[id+1]->apply(sol->thrusts[turn+DEPTH], sol->angles[turn+DEPTH]);
    }

    void move() {
        move(&sol);
    }

    void solve(float time, bool with_seed = false) {
        Solution best;
        if (with_seed) {
            best = sol;
            best.shift();
        } else {
            best.randomize();
            if (r == 0 && pods[id]->dist(cps[1]) > 4000) best.thrusts[0] = 650;
        }
        get_score(&best);

        Solution child;
        while (TIME < time) {
            best.mutate(&child);
            if (get_score(&child) > get_score(&best)) best = child;
        }
        sol = best;
    }

    float get_score(Solution* sol) {
        if (sol->score == -1) {
            vector<float> scores;
            for (Bot* oppBot : oppBots) {
                scores.push_back(get_bot_score(sol, oppBot));
            }

            sol->score = *min_element(scores.begin(), scores.end());
        }

        return sol->score;
    }

    float get_bot_score(Solution* sol, Bot* opp) {
        float score = 0;
        while (turn < DEPTH) {
            move(sol);
            opp->move();
            play();
            if (turn == 0) score += 0.1*evaluate();
            turn++;
        }
        score += 0.9*evaluate();
        load();

        if (r > 0) sols_ct++;

        return score;
    }

    float evaluate() {
        Pod* my_runner = runner(pods[id], pods[id+1]);
        Pod* my_blocker = blocker(pods[id], pods[id+1]);
        Pod* opp_runner = runner(pods[(id+2) % 4], pods[(id+3) % 4]);
        Pod* opp_blocker = blocker(pods[(id+2) % 4], pods[(id+3) % 4]);

        float score = my_runner->score() - opp_runner->score();
        // TODO maybe not a great idea? :)
        score -= my_blocker->dist(my_runner);

        return score;
    }
};

void load() {
    for (int i = 0; i < 4; i++) pods[i]->load();
    turn = 0;
}

void play() {
    float t = 0.0;
    while (t < 1.0) {
        Collision first_col = {NULL, NULL, -1};
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                float col_time = pods[i]->collision_time(pods[j]);
                if (col_time > -1 && col_time + t < 1.0 && (first_col.t == -1 || col_time < first_col.t)) {
                    first_col.a = pods[i];
                    first_col.b = pods[j];
                    first_col.t = col_time;
                }
            }

            // TODO this is wasteful, get rid of it
            float col_time = pods[i]->collision_time(cps[pods[i]->ncpid]);
            if (col_time > -1 && col_time + t < 1.0 && (first_col.t == -1 || col_time < first_col.t)) {
                first_col.a = pods[i];
                first_col.b = cps[pods[i]->ncpid];
                first_col.t = col_time;
            }
        }

        if (first_col.t == -1) {
            for (int i = 0; i < 4; i++) {
                pods[i]->move(1.0 - t);
            }
            t = 1.0;
        } else {
            for (int i = 0; i < 4; i++) {
                pods[i]->move(first_col.t);
            }

            first_col.a->bounce(first_col.b);
            t += first_col.t;
        }
    }

    for (int i = 0; i < 4; i++) {
        pods[i]->end();
    }
}

void print_move(int thrust, float angle, Pod* pod) {
    float a = pod->angle + angle;

    if (a >= 360.0) {
        a = a - 360.0;
    } else if (a < 0.0) {
        a += 360.0;
    }

    a = a * M_PI / 180.0;
    float px = pod->x + cos(a) * 10000.0;
    float py = pod->y + sin(a) * 10000.0;

    char copyright[] = "github.com/inoryy/csb-ai-starter"; // do not remove
    if (thrust == -1) {
        printf("%d %d SHIELD %s\n", (int) round(px), (int) round(py), copyright);
        pod->shield = 4;
    } else if (thrust == 650) {
        pod->has_boost = false;
        printf("%d %d BOOST %s\n", (int) round(px), (int) round(py), copyright);
    } else {
        printf("%d %d %d %s\n", (int) round(px), (int) round(py), thrust, copyright);
    }
}

int main() {
    cin >> laps >> cp_ct;
    for (int i = 0; i < cp_ct; i++) {
        int cx, cy;
        cin >> cx >> cy;
        cps[i] = new Checkpoint(i, cx, cy);
    }

    for (int i = 0; i < 4; i++) pods[i] = new Pod(i);

    pods[0]->partner = pods[1];
    pods[1]->partner = pods[0];
    pods[2]->partner = pods[3];
    pods[3]->partner = pods[2];

    ReflexBot me_reflex;

    SearchBot opp(2);
    opp.oppBots.push_back(&me_reflex);

    SearchBot me;
    me.oppBots.push_back(&opp);

    while (1) {
        r++;

        for (int i = 0; i < 4; i++) {
            int x, y, vx, vy, angle, ncpid;
            cin >> x >> y >> vx >> vy >> angle >> ncpid;
            if (r == 0 && i > 1 && angle > -1) is_p2 = true;
            pods[i]->update(x, y, vx, vy, angle, ncpid);
        }

        now = high_resolution_clock::now();

        float time_limit = r ? 0.142 : 0.98;
        time_limit *= 0.3;

        // use this to test reflex bot behavior
        // me_reflex.move_as_main();

        opp.solve(time_limit*0.15);
        me.solve(time_limit, r > 0);

        if (r > 0) cerr << "Avg iters: " << sols_ct / r << "; Avg sims: " << sols_ct*DEPTH / r << endl;

        print_move(me.sol.thrusts[0], me.sol.angles[0], pods[0]);
        print_move(me.sol.thrusts[DEPTH], me.sol.angles[DEPTH], pods[1]);
    }
}