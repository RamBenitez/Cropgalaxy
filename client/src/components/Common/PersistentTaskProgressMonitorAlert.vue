<script setup lang="ts">
import { library } from "@fortawesome/fontawesome-svg-core";
import { faSpinner } from "@fortawesome/free-solid-svg-icons";
import { FontAwesomeIcon } from "@fortawesome/vue-fontawesome";
import { BAlert, BLink } from "bootstrap-vue";
import { computed, watch } from "vue";

import { TaskMonitor } from "@/composables/genericTaskMonitor";
import { MonitoringRequest, usePersistentProgressTaskMonitor } from "@/composables/persistentProgressMonitor";
import { useShortTermStorage } from "@/composables/shortTermStorage";

library.add(faSpinner);

interface Props {
    monitorRequest: MonitoringRequest;
    useMonitor: TaskMonitor;

    /**
     * The task ID to monitor. Can be a task ID or a short-term storage request ID.
     * If provided, the component will start monitoring the task with the given ID.
     */
    taskId?: string;

    /**
     * If true, the download link will be automatically opened when the task is completed if the user
     * remains on the page.
     *
     * Automatic download will only be possible if the task is completed and the task type is `short_term_storage`.
     */
    enableAutoDownload?: boolean;

    inProgressMessage?: string;
    completedMessage?: string;
    failedMessage?: string;
    requestFailedMessage?: string;
}

const props = withDefaults(defineProps<Props>(), {
    taskId: undefined,
    enableAutoDownload: false,
    inProgressMessage: `Task is in progress. Please wait...`,
    completedMessage: "Task completed successfully.",
    failedMessage: "Task failed.",
    requestFailedMessage: "Request failed.",
});

const { getDownloadObjectUrl } = useShortTermStorage();

const { hasMonitoringData, isRunning, isCompleted, hasFailed, requestHasFailed, storedTaskId, status, start, reset } =
    usePersistentProgressTaskMonitor(props.monitorRequest, props.useMonitor);

const downloadUrl = computed(() => {
    // We can only download the result if the task is completed and the task type is short_term_storage.
    const requestId = props.taskId || storedTaskId;
    if (requestId && props.monitorRequest.taskType === "short_term_storage") {
        return getDownloadObjectUrl(requestId);
    }
    return undefined;
});

if (hasMonitoringData.value) {
    start();
}

watch(
    () => props.taskId,
    (newTaskId, oldTaskId) => {
        if (newTaskId && newTaskId !== oldTaskId) {
            start({ taskId: newTaskId, taskType: props.monitorRequest.taskType, request: props.monitorRequest });
        }
    }
);

watch(
    () => isCompleted.value,
    (completed) => {
        // We check for props.taskId to be defined to avoid auto-downloading when the task is completed and
        // the component is first mounted like when refreshing the page.
        if (completed && props.enableAutoDownload && downloadUrl.value && props.taskId) {
            window.open(downloadUrl.value, "_blank");
        }
    }
);

function dismissAlert() {
    reset();
}
</script>

<template>
    <div v-if="hasMonitoringData" class="d-flex justify-content-end">
        <BAlert v-if="isRunning" variant="info" show>
            <b>{{ inProgressMessage }}</b>
            <FontAwesomeIcon :icon="faSpinner" class="mr-2" spin />
        </BAlert>
        <BAlert v-else-if="isCompleted" variant="success" show dismissible @dismissed="dismissAlert">
            <span>{{ completedMessage }}</span>
            <BLink v-if="downloadUrl" class="download-link" :href="downloadUrl">
                <b>Download here</b>
            </BLink>
        </BAlert>
        <BAlert v-else-if="hasFailed" variant="danger" show dismissible @dismissed="dismissAlert">
            <span>{{ failedMessage }}</span>
            <span v-if="status">
                Reason: <b>{{ status }}</b>
            </span>
        </BAlert>
        <BAlert v-else-if="requestHasFailed" variant="danger" show dismissible @dismissed="dismissAlert">
            <b>{{ requestFailedMessage }}</b>
        </BAlert>
    </div>
</template>
