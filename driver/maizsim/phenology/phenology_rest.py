        # self.order = {
        #     g: [e, li, la, ti, gstt, gddt_ae],
        #     li: [ati],
        #     la: [s],
        #     s: [gf],
        # }

    def update(self, T):
        for stage in queue:
            if stage.ready():
                stage.update(T)
            if stage.over():
                next_stages = self.order[stage]
                self.queue.extend(next_stages)
                self.queue.remove(stage)
                stage.finish()

                #FIXME just for debugging & compatibility
                print("* {}: GDDsum = {}, Growing season T = {}" % (self.gdd_recorder.rate, self.gst_recorder.rate)





    def update(self):
        if self.leaves_appeared < 9:
            T_cur = wthr.soilT
        else:
            T_cur = wthr.airT
        T_cur = max(0, T_cur)

        # converting minute to day decimal, 1= a day
        #dt = initInfo.timeStep/(24*60)

        if not self.germinated:
            self.germinate()
        else:
            T_grow_sum += T_cur
            steps += 1
            # mean growing season temperature since germination, SK 1-19-12
            T_grow = T_grow_sum / steps

            if not self.emerged:
                self.emerge()

            if not self.tassel_initiated:
                self.tassel_initiate()
            else:
                # to be used for C partitoining time scaling, see Plant.cpp
                phyllochrons_from_TI += self.beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil)

            if leaves_appeared < int(leaves_initiated):
                leaves_appeared += self.beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil)
                leaves_appeared = min(leaves_appeared, int(leaves_initiated))

            if self.tassel_initiated and leaves_appeared >= int(leaves_initiated) and not self.silking_done:
                self.anthesis()

            if self.silking_done:
                self.grain_fill()

        d_GTI = self.calc_GTI(T_cur, self.silking_done) * dt
        d_GDD = self.calc_GDD(T_cur) * dt
        GDD_sum += d_GDD
        GTI_sum += d_GTI

        if GDD_sum >= GDD_rating and not self.matured:
            self.mature()


class Development:
    def __init__(self, info):
        leaves_initiated = 0
        leaves_appeared = 0
        leaves_expanded = 0
        anthesis = 0

        germination_rate = 0.
        emergence_rate = 0.

        GDD_sum = 0.
        GDD_grain = 0.
        d_GDD = 0.

        GDD_rating = 1000.
        emerge_gdd = 0.

        P2 = 0.
        phyllochrons_from_TI = 0.

        GTI_sum = 0.
        d_GTI 0.

        Rax_LIR = info.Rmax_LIR
        Rmax_LTAR = info.Rmax_LTAR
        day_length_sensitivity = info.DayLengthSensitive
        Rmax_germination = 0.
        RMax_emergence = 0.

        T_base = 8.0
        T_opt = 30.0
        T_ceil = 40.0

        T_calibration = 20.0

        total_leaf_number = info.genericLeafNo
        juvenile_leaf_number = info.genericLeafNo
        initiated_leaf_number = 5
        youngest_leaf_number = 5
        current_leaf_number = 1

        leaves_at_TI = 1

        init_info = info

        leaves_to_induce = 0
        induction_period = 0.
        inductions = 0

        # converting minute to day decimal, 1 = a day
        dt = info.timeStep / MINUTESPERDAY

        T_grow_sum = 0.
        steps = 0
        T_grow = -99
        T_ind = -99

        phyllochrons_to_silk = info.PhyllochronsToSilk

        # setParms()

        # max rate of germination per day, assume it takes two day at 31 C, needs to be replaced
        Rmax_germination = 0.45*dt
        Rmax_emergence = 0.2388*dt

        # Kim et al. (2007); Kim and Reddy (2004), 0.581 from Yan and Hunt (1999), equivalent phyllochron in CERES
        Rmax_LTAR = Rmax_LTAR*dt

        # cdt changed from 0.524 to test for colorado data used 0.374
        # best fit of K and W (1983), Kim and Reddy (2004) originally 0.978 (for colorado tried 0.558
        Rmax_LIR = Rmax_LIR*dt

        # These Topt and Tceil values from Kim et al. (2007), also see Kim and Reddy (2004), Yan and Hunt (1999), SK
        T_base = 8.0
        T_opt  = 32.1
        T_ceil = 43.7

        leaves_initiated = initiated_leaf_number

        gdd_rating = info.GDD_rating

        P2 = 0.5

    def update(self, weather):
        if self.leaves_appeared < 9:
            T_cur = wthr.soilT
        else:
            T_cur = wthr.airT
        T_cur = max(0, T_cur)

        # converting minute to day decimal, 1= a day
        #dt = initInfo.timeStep/(24*60)

        if not self.germinated:
            #TODO implement germination rate model of temperature.
            # for now assume it germinates immidiately after sowing
            self.germinate()
        else:
            T_grow_sum += T_cur
            steps += 1
            # mean growing season temperature since germination, SK 1-19-12
            T_grow = T_grow_sum / steps

            if not self.emerged:
                self.emerge()

            if not self.tassel_initiated:
                self.tassel_initiate()
            else:
                # to be used for C partitoining time scaling, see Plant.cpp
                phyllochrons_from_TI += self.beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil)

            if leaves_appeared < int(leaves_initiated):
                leaves_appeared += self.beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil)
                leaves_appeared = min(leaves_appeared, int(leaves_initiated))

            if self.tassel_initiated and leaves_appeared >= int(leaves_initiated) and not self.silking_done:
                self.anthesis()

            if self.silking_done:
                self.grain_fill()

        d_GTI = self.calc_GTI(T_cur, self.silking_done) * dt
        d_GDD = self.calc_GDD(T_cur) * dt
        GDD_sum += d_GDD
        GTI_sum += d_GTI

        if GDD_sum >= GDD_rating and not self.matured:
            self.mature()

    def germinate(self):
        #TODO implement germination rate model of temperature.
        # for now assume it germinates immidiately after sowing
        self.germination_rate += self.beta_fn(T_cur, Rmax_Germination, T_opt, T_ceil)
        if (self.germination_rate >= 0.5:
            #TODO record event
            germination.done = true;
            germination.daytime = wthr.daytime;

            print("* Germinated: GDDsum = {}, time step (min) = {}" % (GDD_sum, dt*(24*60)))

    def emerge(self):
        self.emergence_rate += self.beta_fn(T_cur, Rmax_Emergence, T_opt, T_ceil)
        # corn->LvsAppeared = corn->EmergenceRate;
        if self.emergence_rate >= 1.0:
            # corn->LvsAppeared = 1.0

            #TODO record event
            emergence.done = true;
            emergence.daytime = wthr.daytime;

            # reset GDDsum from emergernce, SK
            GDD_sum = 0.0

            # gdd at emergence YY 4/2/09
            emerge_GDD = GDD_sum

            print("* Emergence: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))

    def tassel_initiate(self):
        self.leaves_initiated += self.beta_fn(T_cur, Rmax_LIR, T_opt, T_ceil)
        self.current_leaf_number = (int) self.leaves_initiated

        if self.leaves_initiated >= self.juvenile_leaf_number:
            # inductive phase begins after juvenile stage and ends with tassel initiation

            # Equation 4 in Grant 1989 Ag. J. (81)
            #dt 12/11/2012 broke the equation in two to separate temperature and daylenght effects
            #addLeaf = max(0, 0.1*(juvLeafNo-10.0)*(wthr.dayLength-12.5) + (13.9-1.89*T_cur+0.0795*T_cur*T_cur - 0.001*T_cur*T_cur*T_cur))

            if T_ind == -99:
                # mean temperature during induction period
                T_ind = T_grow

            addLeafTemperature = max(0., 13.6 - 1.89*T_ind + 0.081*T_ind**2 - 0.001*T_ind**3)
            if DayLengthSensitive:
                addLeafPhotoPeriod = max(0., 0.1 * (juvLeafNo - 10.0) * (wthr.dayLength - 12.5))
            else:
                addLeafPhotoPeriod = 0.
            addLeafTotal = addLeafTemperature + addLeafPhotoPeriod

            # effect of photoperiod and temperature on leaf no. used as Grant (1989)
            # Added back the temperature effect on leaf number and revised the algorithm to accumulate addLeafNo to totLeafNo.
            # Changed to respond to mean growing season temperature upto this point.
            # This has little mechanistic basis. Needs improvements. SK 1-19-12
            leaves_to_induce = (leaves_to_induce * inductions + addLeafTotal) / (inductions+1)
            T_ind = (T_ind*inductions + T_cur)/(inductions +1)
            inductions += 1
            inductionPeriod += dt

            # totLeafNo = juvLeafNo + addedLvs/inductionPeriod; //get a mean value over this period
            # LvsAtTI = LvsInitiated; //Should be LvsInitiated. Already confirmed with Soo. 7/27/2006
            # uncomment the following for debugging
            #cout << "* Inductive phase: " << LvsInitiated << " " << totLeafNo << " " << juvLeafNo << " " << addedLvs/inductionPeriod << endl;
            actual_added_leaves = leaves_initiated - juvenile_leaf_number
            if actual_added_leaves >= leaves_to_induce:
                youngestLeaf = totLeafNo = (int) LvsInitiated
                curLeafNo = youngestLeaf

                #TODO record event
                tasselInitiation.done = True
                tasselInitiation.daytime = wthr.daytime
                self.tassel_initiated = True

                leaves_initiated = youngest_leaf
                leaves_at_TI = self.leaves_appeared
                print("* Tassel initiation: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))

    def anthesis(self):
        # Assume 75% Silking occurs at total tip appeared + 3 phyllochrons
        self.anthesis += self.beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil)
        if self.anthesis >= phyllochrons_to_silk: # was 3
            #TODO: record event
            silking.done = True
            silking.daytime = wthr.daytime
            print("* Silking: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))

    def grain_fill(self):
        self.GDD_grain += self.calc_GDD(T_cur) * dt
        if self.GDD_grain >= 170 and not self.begin_grain_fill:
            #TODO GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
            #TODO record event
            beginGrainFill.done = True
            beginGrainFill.daytime = wthr.daytime
            print("* Grain filling begins: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))

    #######
    #######

    def beta_fn(self, T, R_max, To, Tc):
        # beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
        Tb = 0.
        if T <= Tb or Tc <= T:
            return 0.
        if To <= Tb or Tc <= To:
            return 0.
        Tc = Tc*1.
        To = To*1.

        f = (T - Tb) / (To - Tb)
        g = (Tc - T) / (Tc - To)
        beta = 1.0
        alpha = beta * (To - Tb) / (Tc - To)
        return R_max * f**alpha * g**beta


    def calc_GTI(self, T_avg, silked):
        # General Thermal Index, Stewart et al. (1998)
        # Phenological temperature response of maize. Agron. J. 90: 73-79.
        #b1 = 0.0432
        b1 = 0.011178

        #T_opt = 32.2
        #return b1*T_avg**2 * (1 - 0.6667 * T_avg/T_opt)
        #if not Silked:
            #return b1*T_avg**2 * (1 - 0.6667 * T_avg/T_opt)
        #else:
        return b1*T_avg**2 + 5.358

    def calc_GDD(self, T_avg):
        # GDD model with base 8. See Birch et al. (2003) Eu J Agron
        T_base = 8.0
        T_opt = 34.0
        #T_opt = 30.0
        return min(T_avg, T_opt) - T_base
