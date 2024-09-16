#include "elasticRibbonForce.h"

elasticRibbonForce::elasticRibbonForce(elasticRod &m_rod, timeStepper &m_stepper)
{
	rod = &m_rod;
	stepper = &m_stepper;

    gradKappa1 = MatrixXd::Zero(rod->nv,11);
    gradKappa2 = MatrixXd::Zero(rod->nv,11);
    relevantPart = MatrixXd::Zero(11, 2);;
    DDkappa1 = MatrixXd::Zero(11,11);
    DDkappa2 = MatrixXd::Zero(11,11);

    Jbb_2 = MatrixXd::Zero(11,11);
    
    D2kappa1De2.setZero(3,3);
    D2kappa1Df2.setZero(3,3);
    D2kappa1DeDf.setZero(3,3);
    D2kappa2De2.setZero(3,3);
    D2kappa2Df2.setZero(3,3);
    D2kappa2DeDf.setZero(3,3);
    kappa11 = VectorXd::Zero(rod->nv);
    kappa22 = VectorXd::Zero(rod->nv);

    gradTwist = MatrixXd::Zero(rod->nv,11);
    deltam = VectorXd::Zero(rod->ne);
    theta_f = VectorXd::Zero(rod->ne);
    theta_e = VectorXd::Zero(rod->ne);

    DDtwist.setZero(11,11);

    Jtt_2 = MatrixXd::Zero(11,11);

    gradTwistLocal.setZero(11);

    Id3<<1,0,0,
         0,1,0,
         0,0,1;

    fBending_2 = VectorXd::Zero(11);

    fTwisting_2 = VectorXd::Zero(11);

    outputKappa = VectorXd::Zero(rod->nv);
    outputTau = VectorXd::Zero(rod->nv);

    forceVec = VectorXd::Zero(rod->ndof);

    YoungM = rod->youngM;
    ribbonWidth = rod->width;
    ribbonThickness = rod->rodRadius;

    EIMat << 1,0,
             0,1;

    EIMat(0, 0) = 400.0 * YoungM * ribbonWidth * pow(ribbonThickness, 3) / 12;
    EIMat(1, 1) = 1.000 * YoungM * ribbonWidth * pow(ribbonThickness, 3) / 12;

    GJ = 1.4562 * YoungM * ribbonWidth * pow(ribbonThickness, 3) / 12;
}

elasticRibbonForce::~elasticRibbonForce()
{
	;
}

void elasticRibbonForce::computeFribbon()
{
    forceVec = VectorXd::Zero(rod->ndof);

    // bending force
    for(int i=1;i<rod->ne;i++)
    {
        norm_e = rod->edgeLen(i-1);
        norm_f = rod->edgeLen(i);
        te = rod->tangent.row(i-1);
        tf = rod->tangent.row(i);
        d1e = rod->m1.row(i-1);
        d2e = rod->m2.row(i-1);
        d1f = rod->m1.row(i);
        d2f = rod->m2.row(i);

        chi = 1.0 + te.dot(tf);
        tilde_t = (te+tf)/chi;
        tilde_d1 = (d1e+d1f)/chi;
        tilde_d2 = (d2e+d2f)/chi;

        kappa1 = rod->kappa(i,0);
        kappa2 = rod->kappa(i,1);

        Dkappa1De = (1.0/norm_e)*(-kappa1*tilde_t + tf.cross(tilde_d2));
        Dkappa1Df = (1.0/norm_f)*(-kappa1*tilde_t - te.cross(tilde_d2));
        Dkappa2De = (1.0/norm_e)*(-kappa2*tilde_t - tf.cross(tilde_d1));
        Dkappa2Df = (1.0/norm_f)*(-kappa2*tilde_t + te.cross(tilde_d1));

        gradKappa1.row(i).segment(0,3)=-Dkappa1De;
        gradKappa1.row(i).segment(4,3)= Dkappa1De - Dkappa1Df;
        gradKappa1.row(i).segment(8,3)= Dkappa1Df;

        gradKappa2.row(i).segment(0,3)=-Dkappa2De;
        gradKappa2.row(i).segment(4,3)= Dkappa2De - Dkappa2Df;
        gradKappa2.row(i).segment(8,3)= Dkappa2Df;

        kbLocal = (rod->kb).row(i);

        gradKappa1(i,3)=-0.5*kbLocal.dot(d1e);
        gradKappa1(i,7)=-0.5*kbLocal.dot(d1f);
        gradKappa2(i,3)=-0.5*kbLocal.dot(d2e);
        gradKappa2(i,7)=-0.5*kbLocal.dot(d2f);
    }

    for (int i=1; i < rod->ne; i++)
    {
        ci = 4*i-4;
        relevantPart.col(0) = gradKappa1.row(i);
        relevantPart.col(1) = gradKappa2.row(i);
        kappaL = (rod->kappa).row(i) - (rod->kappaBar).row(i);

        fBending_2 = - relevantPart * EIMat * kappaL / rod->voronoiLen(i);

        for (int k = 0; k < 11; k++)
        {
            ind = ci + k;

            forceVec(ind) = forceVec(ind) + fBending_2[k];

            stepper->addForce(ind, -fBending_2[k]); // subtracting elastic force
        }
    }    
    
    // twisting force
    for (int i=0; i < rod->ne; i++)
    {
        theta_f(i) = rod->x(4*i+3);
    }

    for (int i=0; i<rod->ne; i++)
    {
        if (i==0)
        {
            theta_e(i)=0;
        }
        else
        {
            theta_e(i)=theta_f(i-1);
        }
    }

    deltam = theta_f-theta_e;

    for(int i=1; i < rod->ne; i++)
    {
        norm_e = rod->edgeLen(i-1);
        norm_f = rod->edgeLen(i);
        gradTwist.row(i).segment(0,3) = -0.5 / norm_e * (rod->kb).row(i);
        gradTwist.row(i).segment(8,3) = 0.5 / norm_f * (rod->kb).row(i);
        gradTwist.row(i).segment(4,3) = -(gradTwist.row(i).segment(0,3)+gradTwist.row(i).segment(8,3));
        gradTwist(i, 3) = -1;
        gradTwist(i, 7) =  1;
    }

    for(int i=1; i < rod->ne; i++)
    {
        // 2
        value = GJ / rod->voronoiLen(i) * (deltam(i) + rod->refTwist (i) - rod->undeformedTwist(i));
        ci = 4*i-4;
        fTwisting_2 = - value * gradTwist.row(i);
        
        for (int k = 0; k < 11; k++)
        {
            ind = ci + k;
            
            forceVec(ind) = forceVec(ind) + fTwisting_2[k];

            stepper->addForce(ind, - fTwisting_2[k]);
        }
    }

    totalEnergy = 0.0;

    for (int i = 1; i < rod->nv - 1; i++)
    {
        kappaL = (rod->kappa).row(i) - (rod->kappaBar).row(i);

        tau = rod->tau(i);

        totalEnergy = totalEnergy + 0.5 * EIMat(1, 1) * kappaL(1) * kappaL(1) * rod->voronoiLen(i) + 0.5 * GJ * tau * tau * rod->voronoiLen(i);
    }
}

void elasticRibbonForce::computeJribbon()
{
	for(int i=1; i < rod->ne;i++)
    {
        norm_e = rod->edgeLen(i-1);
        norm_f = rod->edgeLen(i);
        te = rod->tangent.row(i-1);
        tf = rod->tangent.row(i);

        d1e = rod->m1.row(i-1);
        d2e = rod->m2.row(i-1);
        d1f = rod->m1.row(i);
        d2f = rod->m2.row(i);

        norm2_e=norm_e*norm_e;
        norm2_f=norm_f*norm_f;

        chi=1.0+te.dot(tf);
        tilde_t=(te+tf)/chi;
        tilde_d1=(d1e+d1f)/chi;
        tilde_d2=(d2e+d2f)/chi;

        kappa1 = rod->kappa(i,0);
        kappa2 = rod->kappa(i,1);

        kbLocal = (rod->kb).row(i);
		
        tt_o_tt = tilde_t * tilde_t.transpose();

        crossMat(tilde_d1,tilde_d1_3d);
        crossMat(tilde_d2,tilde_d2_3d);

        tmp = tf.cross(tilde_d2);
        tf_c_d2t_o_tt=tmp*tilde_t.transpose();
        tt_o_tf_c_d2t=tf_c_d2t_o_tt.transpose();
        kb_o_d2e=kbLocal*d2e.transpose();
        d2e_o_kb=kb_o_d2e.transpose();

        D2kappa1De2=
        1.0/norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t)
        - kappa1 / (chi * norm2_e) * (Id3 - te * te.transpose())
        + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

        tmp=te.cross(tilde_d2);
        te_c_d2t_o_tt = tmp * tilde_t.transpose();
        tt_o_te_c_d2t= te_c_d2t_o_tt.transpose();
        kb_o_d2f= kbLocal * d2f.transpose();
        d2f_o_kb= kb_o_d2f.transpose();

        D2kappa1Df2=
        1.0 / norm2_f * (2* kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t)
        - kappa1 / (chi * norm2_f) * (Id3 - tf*tf.transpose())
        + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

        D2kappa1DeDf=
        - kappa1/(chi* norm_e * norm_f) *(Id3 + te*tf.transpose())
        + 1.0 / (norm_e*norm_f) * (2*kappa1 *tt_o_tt - tf_c_d2t_o_tt +
            tt_o_te_c_d2t - tilde_d2_3d);
        D2kappa1DfDe = D2kappa1DeDf.transpose();

        tmp = tf.cross(tilde_d1);
        tf_c_d1t_o_tt = tmp*tilde_t.transpose();
        tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
        kb_o_d1e = kbLocal * d1e.transpose();
        d1e_o_kb = kb_o_d1e.transpose();

        D2kappa2De2
        = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t)
        - kappa2 / (chi * norm2_e) * (Id3 - te*te.transpose())
        - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

        tmp = te.cross(tilde_d1);
        te_c_d1t_o_tt = tmp*tilde_t.transpose();
        tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
        kb_o_d1f = kbLocal*d1f.transpose();
        d1f_o_kb = kb_o_d1f.transpose();

        D2kappa2Df2
        = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t)
        - kappa2 / (chi * norm2_f) * (Id3 - tf*tf.transpose())
        - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

        D2kappa2DeDf
        = -kappa2/(chi * norm_e * norm_f) * (Id3 + te*tf.transpose())
        + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt 
        - tt_o_te_c_d1t + tilde_d1_3d);

        D2kappa2DfDe = D2kappa2DeDf.transpose();

        D2kappa1Dthetae2 = -0.5 * kbLocal.dot(d2e);
        D2kappa1Dthetaf2 = -0.5 * kbLocal.dot(d2f);
        D2kappa2Dthetae2 =  0.5 * kbLocal.dot(d1e);
        D2kappa2Dthetaf2 =  0.5 * kbLocal.dot(d1f);

        D2kappa1DeDthetae
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d1e)) * tilde_t 
        - 1.0 / chi * (tf.cross(d1e)));
        D2kappa1DeDthetaf
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d1f)) * tilde_t
        - 1.0 / chi * (tf.cross(d1f)));
        D2kappa1DfDthetae
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d1e)) * tilde_t 
        + 1.0 / chi * (te.cross(d1e)));
        D2kappa1DfDthetaf
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d1f)) * tilde_t 
        + 1.0 / chi * (te.cross(d1f)));
        D2kappa2DeDthetae
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d2e)) * tilde_t
        - 1.0 / chi * (tf.cross(d2e)));
        D2kappa2DeDthetaf
        = 1.0 / norm_e * ((0.5 * kbLocal.dot(d2f)) * tilde_t 
        - 1.0 / chi * (tf.cross(d2f)));
        D2kappa2DfDthetae
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d2e)) * tilde_t
        + 1.0 / chi * (te.cross(d2e)));
        D2kappa2DfDthetaf
        = 1.0 / norm_f * ((0.5 * kbLocal.dot(d2f)) * tilde_t
        + 1.0 / chi * (te.cross(d2f)));

        DDkappa1.block(0,0,3,3) =   D2kappa1De2;
        DDkappa1.block(0,4,3,3) = - D2kappa1De2 + D2kappa1DeDf;
        DDkappa1.block(0,8,3,3) =               - D2kappa1DeDf;
        DDkappa1.block(4,0,3,3) = - D2kappa1De2                + D2kappa1DfDe;
        DDkappa1.block(4,4,3,3) =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
        DDkappa1.block(4,8,3,3) =                 D2kappa1DeDf                - D2kappa1Df2;
        DDkappa1.block(8,0,3,3) =                              - D2kappa1DfDe;
        DDkappa1.block(8,4,3,3) =                                D2kappa1DfDe - D2kappa1Df2;
        DDkappa1.block(8,8,3,3) =                                               D2kappa1Df2;

        DDkappa1(3, 3) =   D2kappa1Dthetae2;
        DDkappa1(7, 7) =   D2kappa1Dthetaf2;

        DDkappa1.col(3).segment(0,3) = - D2kappa1DeDthetae;
        DDkappa1.col(3).segment(4,3) =   D2kappa1DeDthetae - D2kappa1DfDthetae;
        DDkappa1.col(3).segment(8,3) =                       D2kappa1DfDthetae;
        DDkappa1.row(3).segment(0,3) =   DDkappa1.col(3).segment(0,3).transpose();
        DDkappa1.row(3).segment(4,3) =   DDkappa1.col(3).segment(4,3).transpose();
        DDkappa1.row(3).segment(8,3) =   DDkappa1.col(3).segment(8,3).transpose();

        DDkappa1.col(7).segment(0,3) = - D2kappa1DeDthetaf;
        DDkappa1.col(7).segment(4,3) =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
        DDkappa1.col(7).segment(8,3) =                       D2kappa1DfDthetaf;
        DDkappa1.row(7).segment(0,3) =   DDkappa1.col(7).segment(0,3).transpose();
        DDkappa1.row(7).segment(4,3) =   DDkappa1.col(7).segment(4,3).transpose();
        DDkappa1.row(7).segment(8,3) =   DDkappa1.col(7).segment(8,3).transpose();

        DDkappa2.block(0,0,3,3) =   D2kappa2De2;
        DDkappa2.block(0,4,3,3) = - D2kappa2De2 + D2kappa2DeDf;
        DDkappa2.block(0,8,3,3) =               - D2kappa2DeDf;
        DDkappa2.block(4,0,3,3) = - D2kappa2De2                + D2kappa2DfDe;
        DDkappa2.block(4,4,3,3) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
        DDkappa2.block(4,8,3,3) =                 D2kappa2DeDf                - D2kappa2Df2;
        DDkappa2.block(8,0,3,3) =                              - D2kappa2DfDe;
        DDkappa2.block(8,4,3,3) =                                D2kappa2DfDe - D2kappa2Df2;
        DDkappa2.block(8,8,3,3) =                                               D2kappa2Df2;

        DDkappa2(3, 3)     = D2kappa2Dthetae2;
        DDkappa2(7, 7)     = D2kappa2Dthetaf2;

        DDkappa2.col(3).segment(0,3) = - D2kappa2DeDthetae;
        DDkappa2.col(3).segment(4,3) =   D2kappa2DeDthetae - D2kappa2DfDthetae;
        DDkappa2.col(3).segment(8,3) =                       D2kappa2DfDthetae;
        DDkappa2.row(3).segment(0,3) =   DDkappa2.col(3).segment(0,3).transpose();
        DDkappa2.row(3).segment(4,3) =   DDkappa2.col(3).segment(4,3).transpose();
        DDkappa2.row(3).segment(8,3) =   DDkappa2.col(3).segment(8,3).transpose();

        DDkappa2.col(7).segment(0,3) = - D2kappa2DeDthetaf;
        DDkappa2.col(7).segment(4,3) =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
        DDkappa2.col(7).segment(8,3) =                       D2kappa2DfDthetaf;
        DDkappa2.row(7).segment(0,3) =   DDkappa2.col(7).segment(0,3).transpose();
        DDkappa2.row(7).segment(4,3) =   DDkappa2.col(7).segment(4,3).transpose();
        DDkappa2.row(7).segment(8,3) =   DDkappa2.col(7).segment(8,3).transpose();

        crossMat(te,teMatrix);

        D2mDe2 = -0.25 / norm2_e * (kbLocal * (te+tilde_t).transpose()
            + (te+tilde_t) * kbLocal.transpose());
        D2mDf2 = -0.25 / norm2_f * (kbLocal * (tf+tilde_t).transpose()
            + (tf+tilde_t) * kbLocal.transpose());
        D2mDeDf = 0.5  / (norm_e*norm_f) * (2.0 / chi * teMatrix
            - kbLocal*tilde_t.transpose());
        D2mDfDe = D2mDeDf.transpose();

        DDtwist.block(0,0,3,3) = D2mDe2;
        DDtwist.block(0,4,3,3) =-D2mDe2 + D2mDeDf;
        DDtwist.block(4,0,3,3) =-D2mDe2 + D2mDfDe;
        DDtwist.block(4,4,3,3) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
        DDtwist.block(0,8,3,3) =-D2mDeDf;
        DDtwist.block(8,0,3,3) =-D2mDfDe;
        DDtwist.block(8,4,3,3) = D2mDfDe - D2mDf2;
        DDtwist.block(4,8,3,3) = D2mDeDf - D2mDf2;
        DDtwist.block(8,8,3,3) = D2mDf2;

        // compute Bending
        len = rod->voronoiLen(i);
        relevantPart.col(0)=gradKappa1.row(i);
        relevantPart.col(1)=gradKappa2.row(i);

        Jbb_2 = MatrixXd::Zero(11,11);

        Jbb_2 = - 1.0/len * relevantPart * EIMat * relevantPart.transpose();

        kappaL = (rod->kappa).row(i) - (rod->kappaBar).row(i);

        temp = - 1.0 / len * kappaL.transpose() * EIMat;
        
        Jbb_2 = Jbb_2 + temp(0) * DDkappa1 + temp(1) * DDkappa2;
        
        // compute Twisting
        Jtt_2 = MatrixXd::Zero(11,11);

        gradTwistLocal = gradTwist.row(i);

        double tw = deltam(i)+rod->refTwist(i) - rod->undeformedTwist(i);

        Jtt_2 = - 1 / len * GJ * (tw * DDtwist + gradTwistLocal * gradTwistLocal.transpose());
        
		for (int j = 0;j < 11;j++)
        {
            for (int k = 0;k < 11;k++)
            {
				ind1 = 4 * i - 4 + j;
				ind2 = 4 * i - 4 + k;
				
                stepper->addJacobian(ind1, ind2, - Jbb_2(k,j) - Jtt_2(k,j) );
            }
        }
    }
}

// Utility
void elasticRibbonForce::crossMat(const Vector3d &a,Matrix3d &b)
{
	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;
}
