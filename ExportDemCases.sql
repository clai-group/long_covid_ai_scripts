with cte_loyalty_score_per_person as (
--isolating the loyalty score here since the sqlserver and oracle versions may have different layouts and table names
select * from loyalty_score_by_patient
where predicted_score > 0.3
),
cte_charlson_per_person as (
--isolating the loyalty score here since the sqlserver and oracle versions may have different layouts and table names
select * from loyalty_cohort_charlson
)
--select count(*) from cte_charlson_per_person; --1246291
--select count(*) from cte_loyalty_score_per_person; --587008 > 0.3 
,
cte_loyal_pasc_cohort as (
-- COVID Positive patients from your ENACT PATIENT_SET_COLLECTION where query name was
-- and LOYALTY SCORE is > 0.3 (?)
select ps.patient_num, c.age, c.sex as sex , c.charlson_index, pd.race_cd as race_cd,oft.concept_cd as ethnicity_cd
from qt_patient_set_collection ps
join cte_charlson_per_person c on c.patient_num = ps.patient_num
join patient_dimension pd on pd.patient_num = ps.patient_num
join observation_fact  oft on oft.patient_num=ps.patient_num
where ps.result_instance_id = 69781 --FILL_IN_PATIENT_SET_RESULT_INSTANCE_ID --9097 who_case, who_control, who pre_control
and oft.concept_cd like 'DEM|HISP:%'
)
select PATIENT_NUM as "patient_num",AGE as "age",sex as "sex_cd",CHARLSON_INDEX,RACE_CD as "race_cd",ETHNICITY_CD as "ethnicity_cd"
--select *
from cte_loyal_pasc_cohort;
