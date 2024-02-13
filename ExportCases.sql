with cte_loyalty_score_per_person as (
--isolating the loyalty score here since the sqlserver and oracle versions may have different layouts and table names
select * from loyalty_score_by_patient
where predicted_score > 0.3
),
cte_charlson_per_person as (
--isolating the loyalty score here since the sqlserver and oracle versions may have different layouts and table names
select * from loyalty_cohort_charlson
)
--select count(*) from cte_charlson_per_person; --2678712
--select count(*) from cte_loyalty_score_per_person; --2678712 > 0.3 1546729
,
cte_loyal_pasc_cohort as (
-- COVID Positive patients from your ENACT PATEINT_SET_COLLECTION where query name was 
-- and LOYALTY SCORE is > 0.3 (?)
select ps.patient_num
from qt_patient_set_collection ps
join cte_loyalty_score_per_person l on l.patient_num = ps.patient_num
join cte_charlson_per_person c on c.patient_num = l.patient_num
where ps.result_instance_id = 9654 --FILL_IN_PATIENT_SET_RESULT_INSTANCE_ID --9097 who_case, who_control, who pre_control
)
--select count(*) from cte_loyal_pasc_cohort; --56705
,
cte_ccsr_map as (
--ccsr map where relevance and phenx are not null
select m.concept_path, cd.concept_cd, ccsr_key, phenx from ccsr_pasc_act_map m
join concept_dimension cd on cd.concept_path = m.concept_path
where relevance is null and phenx is not null
)
--select * from cte_ccsr_map where phenx is null;
,
cte_covid_pos as (
--covid positive codes
select cd.concept_path, cd.concept_cd, 
regexp_replace(substr(cd.concept_cd,instr(cd.concept_cd, ':',1)+1),'\.', '') as ccsr_key,  
'COVID-19' as phenx 
from concept_dimension cd 
where (cd.concept_path like '\ACT\UMLS_C0031437\SNOMED_3947185011\UMLS_C0037088\SNOMED_3947183016\%'
or cd.concept_path like '\ACT\UMLS_C0031437\SNOMED_3947185011\UMLS_C0022885\UMLS_C1335447\%'
or cd.concept_path like '\ACT\UMLS_C0031437\SNOMED_3947185011\UMLS_C0022885\ACT_LOCAL_LAB_ANY_POSITIVE\%')
--where cd.concept_path ='\ACT\UMLS_C0031437\SNOMED_3947185011\UMLS_C0037088\SNOMED_3947183016\ICD10CM_J12.82\' or
--cd.concept_path = '\ACT\UMLS_C0031437\SNOMED_3947185011\UMLS_C0037088\SNOMED_3947183016\ICD10CM_U07.1\'
)
--select * from cte_covid_pos;
,
cte_ccsr_all as (
select concept_path, concept_cd, ccsr_key, phenx from cte_covid_pos
union 
select concept_path, concept_cd, ccsr_key, phenx from cte_ccsr_map)
--select * from cte_ccsr_all where phenx like 'COVID%';
--PATIENT_NUM,START_DATE,CONCEPT_CD,C_FULLNAME,PHENX
select p.patient_num as "patient_num", start_date as "start_date", c.concept_cd as "concept_cd", c.concept_path as "c_fullname", c.phenx as "phenx"  
from  cte_loyal_pasc_cohort p
join observation_fact o on o.patient_num = p.patient_num
join cte_ccsr_all c on c.concept_cd = o.concept_cd
;